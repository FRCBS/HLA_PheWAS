
## -------------------------------------------------------
## Analysis of HLA allele PheWAS results
## -------------------------------------------------------

# libs and functions
source('./src/functions.R')


## -------------------------------------------------------
## read data
## -------------------------------------------------------

## PheWAS results

# discovery
results <- fread('./data/R3_HLA_PheWAS_n_over10_18Feb20.tsv', data.table=F)
results <- results[str_split_fixed(results$locus, ':', 2)[, 2]!='', ] # keep only 4-digit
results <- filter(results, !is.na(p.value))
# replication
results.rep <- fread('./data/R5sansR3_HLA_PheWAS_all_26Mar20.tsv', data.table=F)
results.rep <- results.rep[str_split_fixed(results.rep$locus, ':', 2)[, 2]!='', ] # keep only 4-digit
results.rep <- filter(results.rep, !is.na(p.value))

# phenotype conditional results
results.cond <- fread('./data/R5_HLA_PheWAS_conditional_22Jun20.tsv', data.table=F)
results.cond <- mutate(results.cond, tag=str_split_fixed(pheno, '_', 2)[, 1])
results.cond <- mutate(results.cond, covtag=str_split_fixed(covpheno, '_', 2)[, 1])
results.cond <- results.cond[!is.na(results.cond$p.value), ]
colnames(results.cond)[3] <- 'HLA'

# stats
results.stats           <- fread('./data/R3_allele_stats.tsv', data.table=F)
results.stats$locus     <- formHLA(results.stats$locus)
results.rep.stats       <- fread('./data/R5_allele_stats.tsv', data.table=F)
results.rep.stats$locus <- formHLA(results.rep.stats$locus)


## Previous HLA phewas studies

# Karnes et al.
karnes <- fread('./data/phewas_catalog/hla-phewas-catalog.csv', data.table=F)
karnes <- karnes[nchar(str_split_fixed(karnes$snp, '_', 3)[, 3])==4, ] 
karnes$snp <-  gsub('HLA_', '', karnes$snp, fixed=T)
karnes$snp <-  gsub('_', '*', karnes$snp, fixed=T)
karnes$snp <-  gsub('([0-9]{2})([0-9]{2})$', '\\1:\\2', karnes$snp)

# Liu et al.
liu      <- fread('./data/phewas_catalog/jmedgenet2016.tsv') %>% filter(., grepl('HLA', Gene))
liu$Gene <- gsub('HLA-', '', liu$Gene, fixed=T) 

# Hirate et al.
hirata         <- fread('./data/phewas_catalog/hirata.tsv', data.table=F)
hirata$variant <- gsub('HLA-', '', hirata$variant, fixed=T)


## --------------------------------------------------------
## FDR
## --------------------------------------------------------

results     <- mutate(results, adjusted.p=adaptiveBH(p.value, alpha=0.01, silent=T)$adjPValues) 
results.rep <- mutate(results.rep, adjusted.p=adaptiveBH(p.value, alpha=0.01, silent=T)$adjPValues) 

results.cond     <- mutate(results.cond, adjusted.p=adaptiveBH(p.value, alpha=0.01, silent=T)$adjPValues)
results.cond.fdr <- filter(results.cond, adjusted.p<0.01) 


## --------------------------------------------------------
## Phenotype annotation
## --------------------------------------------------------

# phenotype annotation
fg.anno  <- fread('~/Documents/FinnGen/FINNGEN_ENDPOINTS_DF5_V2_2020-02-11_public.tsv', data.table=F)[, c(1,4,5)]
fg.anno$TAGS <- gsub('#', '', fg.anno$TAGS, fixed=T)
fg.anno$TAGS <- str_split_fixed(fg.anno$TAGS, ',', 10)[, 1]
fg.tag <- fread('~/Documents/FinnGen/FINNGEN_endpoint_file_format_V1.2_Public.tsv', data.table=F) # tag definitions
fg.tag <- unite(fg.tag, DESC, CHAPTER, OTHER, sep='')
fg.tag$DESC <- str_split_fixed(fg.tag$DESC, ' \\(', 2)[, 1]
fg.tag$code <- gsub('#', '', fg.tag$code, fixed=T)

# Join annotation with results
results.anno <- left_join(results, fg.anno, by=c('pheno'='NAME'))
results.anno <- results.anno[, c(7, 1, 8, 2:6)]
colnames(results.anno)[c(1:3)] <- c('tag', 'pheno', 'longname') 
results.anno[is.na(results.anno)] <- ' '
# include replication results
results.rep2d <- left_join(unite(results.anno, phenoallele, pheno, locus, sep='_')[, 1:2],
                           unite(results.rep, phenoallele, pheno, locus, sep='_'), by='phenoallele')[, 3:5]
colnames(results.rep2d) <- c('p.value.replication', 'beta.replication', 'SEbeta.replication')
results.anno <- cbind(results.anno, results.rep2d)
mode(results.anno$beta) <- 'numeric'
mode(results.anno$SEbeta) <- 'numeric'

# join stats
results.stats      <- unite(results.stats, phenoallele, pheno, locus, sep='_')
colnames(results.rep.stats)  <- paste0(colnames(results.rep.stats), '.replication')
results.rep.stats  <- unite(results.rep.stats, phenoallele, pheno.replication, locus.replication, sep='_')
results.anno.stats <- unite(results.anno, phenoallele, pheno, locus, sep='_', remove=F)
results.anno.stats <- left_join(results.anno.stats, results.stats, by='phenoallele')
results.anno.stats <- left_join(results.anno.stats, results.rep.stats, by='phenoallele')
results.anno.stats <- results.anno.stats[, c('tag', 'pheno', 'longname', 'locus', 'p.value', 'beta', 'SEbeta', 'adjusted.p', 
                                             'n', 'cases', 'locus.cases.n', 'imp.pprob.mean', 'imp.pprob.sd', 'p.value.replication', 
                                             'beta.replication', 'SEbeta.replication', 'n.replication', 'cases.replication', 
                                             'locus.cases.n.replication', 'imp.pprob.mean.replication', 'imp.pprob.sd.replication')]
results.anno.stats.fdr <- filter(results.anno.stats, adjusted.p<0.01)


## --------------------------------------------------------
## Filter association results, keep replicated
## --------------------------------------------------------

# select replicated and consistent effect direction
tmp <- results.anno.stats.fdr
tmp$REP <- ifelse( (tmp$beta>0 & !is.na(tmp$beta) & tmp$beta.replication>0 & !is.na(tmp$beta.replication)) | 
                     (tmp$beta<0 & !is.na(tmp$beta) & tmp$beta.replication<0 & !is.na(tmp$beta.replication)), T, F)

results.anno.stats.fdr.replicated <- filter(tmp, p.value.replication<0.01 & REP==T) %>% dplyr::select(., -REP) %>% na.omit


## --------------------------------------------------------
## PheWAS summary
## --------------------------------------------------------

# data matrix for plotting
results.anno.2 <- left_join(results.anno, fg.tag, by=c('tag'='code'))
results.anno.2 <- mutate(results.anno.2, gene=str_split_fixed(locus, '\\*', 2)[, 1])
results.anno.2 <- mutate(results.anno.2, logP=log10(p.value)*(-1))
results.anno.2 <- mutate(results.anno.2, logP2=log10(logP+1))

# edit some pheno names
results.anno.2$DESC[is.na(results.anno.2$DESC)] <- 'Miscellaneous, not yet classified endpoints'
results.anno.2$DESC[results.anno.2$DESC=='Psychiatric endpoints from Katri Räikkönen'] <- 'Psychiatric endpoints'
results.anno.2$DESC[
  results.anno.2$DESC=='XVIII Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified'] <- 
  'XVIII Symptoms, signs and abnormal clinical and laboratory findings'
results.anno.2$DESC[
  results.anno.2$DESC=='III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism'] <- 
  'III Diseases of the blood, blood-forming organs and immune system'

# palette
my.pal <- c("#ce6465", "#5cc247", "#5d3ecd",  "#aab82d","#b64bde","#498433","#db48bc",
            "#58bd90","#522d8e","#d4a13c","#6b72d8","#e34b29","#52b0c2","#cf3341","#407e68","#d74681",
            "#9baa4d","#9d4298","#9bad7a","#45386e","#c86c2e","#77a1da","#822c29","#c18bd2","#375024",
            "#cc899f","#85722f","#516894","#d39670","#773359","#6d4529")

# arrange phenos by overll p-value
tag.beta.sum <- tapply(1:nrow(results.anno.2), results.anno.2$DESC, function(x) {
  tmp <- results.anno.2[x, 'p.value'] 
  (tmp %>% as.numeric %>% abs %>% sum(1/., na.rm=T)) / nrow(results.anno.2[x, ])
  })
tag.beta.sum   <- data.frame(DESC=names(tag.beta.sum), SUM=tag.beta.sum)
results.anno.2 <- left_join(results.anno.2, tag.beta.sum, by='DESC')
results.anno.2 <- arrange(results.anno.2, SUM)
results.anno.2$DESC <- factor(results.anno.2$DESC, levels=results.anno.2$DESC %>% unique)

# manual edit of colors
results.anno.2$col[results.anno.2$DESC=='IV Endocrine, nutritional and metabolic diseases'] <- '#52b0c2'
results.anno.2$col[results.anno.2$DESC=='Diabetes endpoints'] <- '#ff363c'
results.anno.2$col[results.anno.2$DESC=='Rheuma endpoints'] <- '#58bd90' #'#407e68'
results.anno.2[, c('DESC', 'col')] %>% distinct
results.anno.2$REPL <- ifelse(results.anno.2$p.value.replication<0.01, '23', '21')

# subset data for plotting
results.anno.2.1 <- results.anno.2[-(1:which(grepl('mastoid', results.anno.2$DESC))[1]-1), ]
# different pch for signif. thresholds
logp2size <- function(x) ifelse(x<0.919132, 0.1, 0.15)

# Plot top phenotypes
p.hla.phewas <- ggplot(results.anno.2.1, aes(logP2, DESC)) +
  geom_jitter(color=results.anno.2.1$col, fill=results.anno.2.1$col,
              size=logp2size(results.anno.2.1$logP2), alpha=0.5, shape=21, height=0.4) +
  scale_x_continuous(limits=c(0, 2.7), breaks=c(0, 1, 2), expand=c(0, 0)) +
  xlab(expression('log'[10]*'(-log'[10]*'('*italic(p)*')+1)')) + ylab('phenotype category') +
  geom_vline(xintercept=log10((log10(5e-8)*(-1))+1), linetype='dashed', size=0.3, alpha=0.8) +
  geom_vline(xintercept=log10((log10(filter(results.anno, adjusted.p<0.01)$p.value %>% max)*(-1))+1), 
             linetype='dashed', size=0.3, alpha=0.8) +
  facet_wrap(~ gene, nrow=1, ncol=7) +
  theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
        axis.text.x=element_text(size=7, angle=0, hjust=0.5), axis.text.y=element_text(size=8.6), 
        panel.border=element_rect(fill=NA, size=.3),
        axis.title.y=element_text(color='black'), legend.position='none',
        strip.background=element_rect(fill=NA)) 



## --------------------------------------------------------
## comparison with other HLA phewas studies
## --------------------------------------------------------

# Karnes et al.
karnes.comp <- map2(
  c('^type 1 diabetes$', '^type 2 diabetes$', 'celiac disease', '^rheumatoid arthritis$', 'psoriasis', 'ketoacidosis', 
    'intestinal infection', 'other infectious and parasitic diseases', 'pneumonia', 'systemic lupus erythematosus', 
    'systemic sclerosis', 'juvenile rheumatoid arthritis'),
  c('T1D', 'T2D', 'K11_COELIAC', 'M13_RHEUMA', 'psoriasis', 'ketoacidosis', 'AB1_INTESTINAL_INFECTIONS', 'AB1_INFECTIONS',
    'pneumonia', 'M13_SLE', 'M13_SYSTSLCE', 'M13_JUVERHEU|M13_JUVEARTH'),
  function(x, y) {
    tmp1 <- filter(karnes, grepl(x, phewas_string, ignore.case=T), p<1e-3)[, c('snp', 'odds_ratio', 'p')]
    tmp2 <- filter(results.anno.replicated, grepl(y, pheno, ignore.case=T))[, c('locus', 'beta', 'p.value')] 
    tmp2$beta <- exp(tmp2$beta)
    colnames(tmp2) <- colnames(tmp1)
    tmp <- inner_join(tmp1, tmp2, by='snp') %>% na.omit
    colnames(tmp) <- c('snp', 'odds_ratio_study', 'p_study', 'odds_ratio_y', 'p_y')
    if(nrow(tmp)>0) { data.frame(tmp, phenotype=x) } else NA
  })
karnes.comp <- karnes.comp[!is.na(karnes.comp)] %>% do.call(rbind, .) %>% distinct
karnes.comp$phenotype <- gsub('\\$|\\^', '', karnes.comp$phenotype)
karnes.comp$study <- 'Karnes et al. (2017)'

# Liu et al.
liu.comp <- map2(
  c('diabetes with other specified manifestations', 'celiac disease', 'rheumatoid arthritis', 'psoriasis$', 'ketoacidosis', 
    'pneumonia', 'systemic sclerosis', 'lichen planus'),
  c('T1D', 'K11_COELIAC', 'M13_RHEUMA', 'psoriasis', 'ketoacidosis',
    'pneumonia', 'M13_SYSTSLCE', 'L12_LICHENPLANUS'),
  function(x, y) {
    tmp1 <- filter(liu, grepl(x, Phenotype, ignore.case=T), p<1e-3)[, c('Gene', 'OR', 'p')]
    tmp2 <- filter(results.anno.replicated, grepl(y, pheno, ignore.case=T))[, c('locus', 'beta', 'p.value')] 
    tmp2$beta <- exp(tmp2$beta)
    tmp2$locus <- str_split_fixed(tmp2$locus, '\\*', 2)[, 1]
    colnames(tmp2) <- colnames(tmp1)
    tmp2 <- lapply(tmp2$Gene %>% unique, function(x) tmp2[tmp2$Gene==x, ] %>% filter(p==min(p))) %>% do.call(rbind, .)
    tmp <- inner_join(tmp1, tmp2, by='Gene') %>% na.omit %>% distinct
    colnames(tmp) <- c('snp', 'odds_ratio_study', 'p_study', 'odds_ratio_y', 'p_y')
    if(nrow(tmp)>0) { data.frame(tmp, phenotype=x) } else NA
  })
liu.comp <- liu.comp[!is.na(liu.comp)] %>% do.call(rbind, .)
liu.comp$phenotype <- gsub('\\$|\\^', '', liu.comp$phenotype)
liu.comp$phenotype <- gsub("diabetes with other specified manifestations", 'diabetes', liu.comp$phenotype)
liu.comp$study <- 'Liu et al. (2016)'

# Hirate et al.
hirata.comp <- map2(
  c('T1D', 'nephrotic syndrome', 'rheumatoid arthritis', 'asthma'),
  c('T1D', 'N14_NEPHROTICSYND', 'M13_RHEUMA', 'J10_ASTHMA'),
  function(x, y) {
    tmp1 <- filter(hirata, grepl(x, pheno, ignore.case=T), pval<1e-3)[, c('gene', 'Effect', 'pval')]
    tmp2 <- filter(results.anno.replicated, grepl(y, pheno, ignore.case=T))[, c('locus', 'beta', 'p.value')] 
    tmp2$beta <- exp(tmp2$beta)
    tmp2$locus <- str_split_fixed(tmp2$locus, '\\*', 2)[, 1]
    colnames(tmp2) <- colnames(tmp1)
    tmp2 <- lapply(tmp2$gene %>% unique, function(x) tmp2[tmp2$gene==x, ] %>% filter(pval==min(pval))) %>% do.call(rbind, .)
    tmp <- inner_join(tmp1, tmp2, by='gene') %>% na.omit %>% distinct
    colnames(tmp) <- c('snp', 'odds_ratio_study', 'p_study', 'odds_ratio_y', 'p_y')
    if(nrow(tmp)>0) { data.frame(tmp, phenotype=x) } else NA
  })
hirata.comp <- hirata.comp[!is.na(hirata.comp)] %>% do.call(rbind, .)
hirata.comp$phenotype <- gsub('\\$|\\^', '', hirata.comp$phenotype)
hirata.comp$phenotype <- gsub("T1D", 'diabetes', hirata.comp$phenotype)
hirata.comp$study <- 'Hirata et al. (2019)'

# Collect together
studies.comp <- rbind(karnes.comp, liu.comp, hirata.comp)
studies.comp$phenotype <- gsub('type 1 ', '', studies.comp$phenotype)
studies.comp <- unite(studies.comp, group, phenotype, study, sep='_', remove=F)




## --------------------------------------------------------
## shared alleles by phenotype tag
## --------------------------------------------------------

# HLA alleles associated with more than 1 tag
results.alleles <- map(unique(results.anno.stats.fdr.replicated$tag), function(x) { 
  out <- filter(results.anno.stats.fdr.replicated, tag==x)$locus %>% unique
  data.frame(tag=x, HLA=out, stringsAsFactors=F) # uniq HLAs per category
}) %>% do.call(rbind, .) 
results.alleles <- results.alleles %>% add_count(HLA, sort=F, name='hla_n') 
results.alleles <- results.alleles %>% add_count(tag, sort=F, name='tag_n')
results.alleles <- filter(results.alleles, hla_n>1)
write.table(results.alleles, './data/shared.hlas.tsv', sep='\t', quote=F, row.names=F)

# conditional results; p.value filtering
results.cond.f <- filter(results.cond, p.value<5e-8)
results.cond.f$HLA <- formHLA(results.cond.f$HLA)
# write out
write.table(results.cond.f, './results/tables/phenotype_conditional_23Jun20.tsv', quote=F, sep='\t', row.names=F)


# for a given tag, find independently associated HLAs
results.cond.f.tags <- results.cond.f$tag %>% unique %>% map(function(x) {
  tmp <- filter(results.cond.f, tag==x)$HLA %>% table
  tmp <- data.frame(tag=x, tmp)
  colnames(tmp)[2:3] <- c('HLA', 'indep.count')
  tmp
  }) %>% do.call(rbind, .)
results.cond.f.tags <- results.cond.f.tags %>% add_count(HLA, sort=F, name='hla_n') # number of HLA in all tags
results.cond.f.tags <- results.cond.f.tags %>% distinct(., HLA, hla_n, .keep_all=T) %>% filter(hla_n>1)
results.cond.f.tags <- results.cond.f.tags %>% add_count(tag, sort=F, name='tag_n') %>% 
  arrange(hla_n)

# list tags per HLA
results.cond.f.tags.2 <-  map(results.cond.f.tags$HLA %>% unique, function(x) {
  data.frame(HLA=x, tag=filter(results.cond.f, HLA==x)$tag %>% unique)
}) %>% do.call(rbind, .)
results.cond.f.tags.2     <- left_join(results.cond.f.tags.2, results.cond.f.tags %>% dplyr::select(-tag), by='HLA')
results.cond.f.tags.2     <- left_join(results.cond.f.tags.2, fg.tag, by=c('tag'='code'))
results.cond.f.tags.2$HLA <- factor(results.cond.f.tags.2$HLA, 
                                    levels=results.cond.f.tags.2$HLA %>% table %>% sort(decreasing=T) %>% names)
results.cond.f.tags.2$tag <- factor(results.cond.f.tags.2$tag, 
                                    levels=results.cond.f.tags.2$tag %>% table %>% sort(decreasing=T) %>% names)
results.cond.f.tags.2$DESC <- factor(results.cond.f.tags.2$DESC, 
                                     levels=results.cond.f.tags.2$DESC %>% table %>% sort(decreasing=F) %>% names)


# select infectious/autoimmune phenos for each significant HLA
results.cond.f2 <- lapply(results.cond.f$HLA %>% unique, function(i) {
  out <- filter(results.cond.f, HLA==i)
  out <- out[
    grepl('AB1|INFECT|ITIS|SEPSIS', out$pheno) | 
      grepl('G6_DEMYEL|E4_DIABETES|J10_ASTHMA|L12_PSORIASIS|M13_ANKYLOSPON_STRICT|M13_RHEUMA', out$pheno), 
  ]
  out[out$covpheno %in% out$pheno, ]
}) %>% do.call(rbind, .)
# select phenos: infectious independently of autoimmune
results.cond.f2 <- filter(results.cond.f2, grepl('AB1|INFECT|ITIS|SEPSIS', pheno),
                          grepl('G6_DEMYEL|E4_DIABETES|J10_ASTHMA|L12_PSORIASIS|M13_ANKYLOSPON_STRICT|M13_RHEUMA', 
                                covpheno))
# add peno annotation
results.cond.f2 <- left_join(results.cond.f2, results.anno[, 2:3], by='pheno')
results.cond.f2$haplo <- rep(NA, nrow(results.cond.f2))

# arrange HLA alleles by haplotype
results.cond.f2$HLA <- factor(results.cond.f2$HLA, levels=
         c("C*07:01", "B*08:01", "DRB1*03:01", "DQB1*02:01", 
           "B*13:02", "DRB1*04:01", "DQA1*03:01", "DQB1*03:02",
           "DQA1*05:01", "DRB1*13:01", "DQA1*01:03", "DQB1*06:03"))




