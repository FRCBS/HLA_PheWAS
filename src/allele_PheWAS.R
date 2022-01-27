
## -------------------------------------------------------
## Analysis and plotting of HLA allele PheWAS results
## -------------------------------------------------------

# libs and funcs
source('./src/functions.R')


## -------------------------------------------------------
## data
## -------------------------------------------------------

# HLA PheWAS results

# discovery cohort
results <- fread('./data/R3_HLA_PheWAS_n_over10_18Feb20.tsv', data.table=F)
results <- results[str_split_fixed(results$locus, ':', 2)[, 2]!='', ] # keep only 4-digit
results <- filter(results, !is.na(p.value))
# replication cohort
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
r3.freqs                <- fread('data/R3_HLA_allelefreqs.tsv', data.table=F)
r5.freqs                <- fread('data/R5_HLA_allelefreqs.tsv', data.table=F)

# other hla phewas studies
# karnes
karnes <- fread('./data/phewas_catalog/hla-phewas-catalog.csv', data.table=F)
karnes <- karnes[nchar(str_split_fixed(karnes$snp, '_', 3)[, 3])==4, ] 
karnes$snp <-  gsub('HLA_', '', karnes$snp, fixed=T)
karnes$snp <-  gsub('_', '*', karnes$snp, fixed=T)
karnes$snp <-  gsub('([0-9]{2})([0-9]{2})$', '\\1:\\2', karnes$snp)
# liu
liu      <- fread('./data/phewas_catalog/jmedgenet2016.tsv') %>% filter(., grepl('HLA', Gene))
liu$Gene <- gsub('HLA-', '', liu$Gene, fixed=T) 
# hirata
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
## Pheno annotation
## --------------------------------------------------------

# phenotype annotation
fg.anno  <- fread('../..//Documents/FinnGen/FINNGEN_ENDPOINTS_DF5_V2_2020-02-11_public.tsv', data.table=F)[, c(1,4,5, 13)]
fg.anno$TAGS <- gsub('#', '', fg.anno$TAGS, fixed=T)
fg.anno$TAGS <- str_split_fixed(fg.anno$TAGS, ',', 10)[, 1]
fg.tag <- fread('../../Documents/FinnGen/FINNGEN_endpoint_file_format_V1.2_Public.tsv', data.table=F) # tag definitions
fg.tag <- unite(fg.tag, DESC, CHAPTER, OTHER, sep='')
fg.tag$DESC <- str_split_fixed(fg.tag$DESC, ' \\(', 2)[, 1]
fg.tag$code <- gsub('#', '', fg.tag$code, fixed=T)

# ICD
icd.main <- fg.anno[4:18, 2:3]
fg.anno <- data.frame(fg.anno, 
                      map(fg.anno$HD_ICD_10, function(x) {
                        out <- explain_table(x)
                        ICD_category <- filter(icd.main, LONGNAME==out$chapter)$NAME
                        ICD_category <- ifelse(length(ICD_category)==0, NA, ICD_category)
                        if(grepl('^E', x)) ICD_category <- icd.main$NAME[4]
                        if(grepl('^F', x)) ICD_category <- icd.main$NAME[5]
                        out$ICD_category <- ICD_category
                        out %>% return
                      }) %>% do.call(rbind, .)
                      )


# Join annotation with results
results.anno <- left_join(results, fg.anno, by=c('pheno'='NAME'))
results.anno <- results.anno[, c(7, 1, 8, 2:21)]
colnames(results.anno)[c(1:3)] <- c('tag', 'pheno', 'longname') 

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

# replication
tmp <- results.anno.stats.fdr
tmp$REP <- ifelse( (tmp$beta>0 & !is.na(tmp$beta) & tmp$beta.replication>0 & !is.na(tmp$beta.replication)) | 
                     (tmp$beta<0 & !is.na(tmp$beta) & tmp$beta.replication<0 & !is.na(tmp$beta.replication)), T, F)

results.anno.stats.fdr.replicated <- filter(tmp, p.value.replication<0.01 & REP==T) %>% dplyr::select(., -REP) %>% na.omit
results.anno.stats.fdr.replicated$FisherP <- map(1:nrow(results.anno.stats.fdr.replicated), function(x) {
  tryCatch(sumlog(results.anno.stats.fdr.replicated[x, c("p.value", "p.value.replication")])$p, 
           error=function(e) NA) %>% return
})

# HLA allele frequencies in discovery and replication cohorts
results.anno.stats.fdr.replicated <- inner_join(results.anno.stats.fdr.replicated, r3.freqs, 
                                                by=c('locus'='HLA')) %>% inner_join(., r5.freqs, by=c('locus'='HLA'))
results.anno.stats.fdr.replicated <- results.anno.stats.fdr.replicated[, c(1:4, 23, 5:13, 24, 14:22)]
colnames(results.anno.stats.fdr.replicated)[5]  <- 'HLA.allele.freq'
colnames(results.anno.stats.fdr.replicated)[15] <- 'replication.HLA.allele.freq'

# write out
results.anno.stats.list <- tapply(1:nrow(results.anno.stats), results.anno.stats$pheno, function(x) 
  results.anno.stats[x, ]) # to list by pheno

results.anno.stats.fdr.list <- tapply(1:nrow(results.anno.stats.fdr), results.anno.stats.fdr$pheno, 
                                      function(x) results.anno.stats.fdr[x, ]) 

results.anno.stats.fdr.replicated.list <- tapply(1:nrow(results.anno.stats.fdr.replicated), 
                                                 results.anno.stats.fdr.replicated$pheno, function(x) 
                                                   results.anno.stats.fdr.replicated[x, ]) 


write.table(lapply(results.anno.stats.list, function(x) rbind(data.frame(x), rep(' ', ncol(results.anno.stats)))) %>% 
              do.call(rbind, .), './results/R3_PheWAS_full.tsv', sep='\t', row.names=F)

write.table(lapply(results.anno.stats.fdr.list, function(x) rbind(data.frame(x), rep(' ', ncol(results.anno.stats.fdr)))) %>% 
              do.call(rbind, .), './results/R3_PheWAS_fdr001.tsv', sep='\t', row.names=F)

fwrite(map(results.anno.stats.fdr.replicated.list, function(x) 
  rbind(data.frame(x), rep(' ', ncol(results.anno.stats.fdr.replicated)))) %>% do.call(rbind, .) %>% data.frame,
  './results/tables/PheWAS_fdr001_replication.tsv', sep='\t', row.names=F)


## Fisher's method p
results.anno$FisherP <- map(1:nrow(results.anno), function(x) {
  tryCatch(sumlog(results.anno[x, c("p.value", "p.value.replication")])$p, 
           error=function(e) NA) %>% return
})


## --------------------------------------------------------
## stats
## --------------------------------------------------------

# alleles
results$locus %>% unique %>% length
filter(results.anno.stats.fdr.replicated, imp.pprob.mean>0.5, imp.pprob.mean.replication>0.5) %>% .$locus %>% unique %>% length

# number of assocs
results.anno.stats.fdr.replicated %>% nrow
# number of phenotypes
results.anno.stats.fdr.replicated$pheno %>% unique %>% length

# specific alleles and phenos
filter(results.anno.stats.fdr.replicated, grepl('DIAB|T1D|T2D|^D1_|^D2_|^DM|E4_D1_|E4_D2_|E4_DM|INSULIN', pheno)) %>% nrow
filter(results.anno.stats.fdr.replicated, grepl('DRB1', locus)) %>% nrow
filter(results.anno.stats.fdr.replicated, grepl('DQA1', locus)) %>% nrow
filter(results.anno.stats.fdr.replicated, grepl('DQB1', locus)) %>% nrow
filter(results.anno.stats.fdr.replicated, grepl('^B', locus)) %>% nrow
filter(results.anno.stats.fdr.replicated, grepl('^C', locus)) %>% nrow
filter(results.anno.stats.fdr.replicated, grepl('^A', locus)) %>% nrow
filter(results.anno.stats.fdr.replicated, grepl('DPB1', locus)) %>% nrow

# number of HLAs per phenotype
results.anno.stats.fdr.replicated$pheno %>% unique %>% lapply(., function(x) {
  data.frame(pheno=x,
             hla.num=filter(results.anno.stats.fdr.replicated, pheno==x)$locus %>% unique %>% length)
}) %>% do.call(rbind, .) %>% arrange(hla.num) %>% tail



## --------------------------------------------------------
## PheWAS plot
## --------------------------------------------------------

my.pal <- c(my.pal, awtools::mpalette, 
            ghibli_palette("PonyoMedium") %>% as.vector, 
            ghibli_palette("MononokeMedium") %>% as.vector)[1:47]

results.anno$FisherP <- as.numeric(results.anno$FisherP)

tags <- fread('./data/FinnGen_tags.txt', data.table=F)
tags$code <- gsub('#', '', tags$code, fixed=T)
tags$CHAPTER <- tags$CHAPTER %>% str_split_fixed(., ' \\(', 2) %>% .[, 1] 
tags <- tags %>% filter(., CHAPTER!='', code!='ICDMAIN')
tags$CHAPTER[
  tags$CHAPTER=="III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism"] <- 
  "III Diseases of the blood and blood-forming organs"

results.anno.3 <- mutate(results.anno,   logP=log10(FisherP)*(-1))
results.anno.3 <- mutate(results.anno.3, logP2=log10(logP+1))
results.anno.3 <- mutate(results.anno.3, gene=str_split_fixed(locus, '\\*', 2)[, 1])

results.anno.3$tag <- factor(results.anno.3$tag, levels=tags$code)
results.anno.3$ICD_chapter <- tags$CHAPTER[results.anno.3$tag] %>% 
  factor(., levels=tags$CHAPTER)
na.ind <- is.na(results.anno.3$ICD_chapter) 

results.anno.3.na <- results.anno.3[na.ind==T, ]
results.anno.3 <- results.anno.3[na.ind==F, ]

results.anno.3$col <- my.pal[1:length(unique(results.anno.3$ICD_chapter))][results.anno.3$ICD_chapter]

# rename factor levels to include number of phenos in each
icd.chap.new <- map(levels(results.anno.3$ICD_chapter), function(x) {
  out <- filter(results.anno.3, ICD_chapter==x)
  n <- length(unique(out$pheno))
  paste0(x, ' (', n, ')') %>% return
}) %>% unlist
levels(results.anno.3$ICD_chapter) <- icd.chap.new
# reverse order
results.anno.3$ICD_chapter <- fct_rev(results.anno.3$ICD_chapter)

# factorise HLA gene variable
results.anno.3$gene <- factor(results.anno.3$gene, levels=c('A' ,'C', 'B', 'DRB1', 'DQA1', 'DQB1', 'DPB1'))

# plot
logp2size <- function(x) ifelse(x<0.919132, 0.01, 0.1)
p.hla.phewas <- ggplot(results.anno.3, aes(logP2, ICD_chapter)) +
  geom_jitter(color=results.anno.3$col, fill=results.anno.3$col,
              size=logp2size(results.anno.3$logP2), alpha=0.8, shape=19, height=0.4) +
  scale_x_continuous(limits=c(0, 2.7), breaks=c(0, 1, 2), expand=c(0, 0)) +
  xlab(expression('log'[10]*'(-log'[10]*'('*italic(p)*')+1)')) + ylab('ICD chapter') +
  geom_vline(xintercept=log10((log10(5e-8)*(-1))+1), linetype='dashed', size=0.3, alpha=0.8) +
  geom_vline(xintercept=log10((log10(filter(results.anno, adjusted.p<0.01)$p.value %>% max)*(-1))+1), 
             linetype='dashed', size=0.3, alpha=0.8) +
  facet_wrap(~ gene, nrow=1, ncol=7) +
  theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
        axis.text.x=element_text(size=8, angle=0, hjust=0.5, color='black'), 
        axis.text.y=element_text(size=8.6, color='black'), 
        panel.border=element_rect(fill=NA, size=.3),
        axis.title.y=element_text(color='black', size=12), 
        axis.title.x=element_text(color='black', size=11), 
        legend.position='none', strip.background=element_rect(fill=NA),
        strip.text=element_text(size=11)) 



## --------------------------------------------------------
## PheWAS QQ plot
## --------------------------------------------------------

jpeg('./results/figures/FigureS1_qqplot_by_ICD.jpg', height=18.5, width=12, units='in', res=600)  
op <- par(mfrow=c(6,4), oma=c(1, 1, 3, 1))
icd.qqplots <- map(results.anno.3$ICD_chapter %>% levels %>% rev %>% .[1:22], function(x) {
  print(x)
  out <- filter(results.anno.3, ICD_chapter==x)$FisherP %>% na.omit
  qqPlot(out, main=str_wrap(x, 40), cex.main=1)
  mtext("Q-Q plots by ICD chapter", side=3, line=0, adj=0, outer=T, cex=1.1, font=2)
})
dev.off()
par(op)

jpeg('./results/figures/FigureS1_qqplot_by_Gene.jpg', height=7, width=12, units='in', res=600)  
op <- par(mfrow=c(2,4), oma=c(1, 1, 3, 1))
icd.qqplots <- map(results.anno.3$gene %>% levels, function(x) {
  print(x)
  out <- filter(results.anno.3, gene==x)$FisherP %>% na.omit
  qqPlot(out, main=paste0('HLA-', x), cex.main=1.5)
  mtext("Q-Q plots by HLA gene", side=3, line=0, adj=0, outer=T, cex=1.1, font=2)
})
dev.off()
par(op)

# append into one jpg and convert to PDF
system("convert -append ./results/figures/FigureS1_qqplot_*.jpg ./results/figures/FigureS1_qqplots.jpg")
system("img2pdf ./results/figures/FigureS1_qqplots.jpg --output ./results/figures/FigureS1_qqplots.pdf")


## --------------------------------------------------------
## Comparison with other HLA phewas studies
## --------------------------------------------------------

karnes.comp <- map2(
  c('^type 1 diabetes$', '^type 2 diabetes$', 'celiac disease', '^rheumatoid arthritis$', 'psoriasis', 'ketoacidosis', 
    'intestinal infection', 'other infectious and parasitic diseases', 'pneumonia', 'systemic lupus erythematosus', 
    'systemic sclerosis', 'juvenile rheumatoid arthritis'),
  c('T1D', 'T2D', 'K11_COELIAC', 'M13_RHEUMA', 'psoriasis', 'ketoacidosis', 'AB1_INTESTINAL_INFECTIONS', 'AB1_INFECTIONS',
    'pneumonia', 'M13_SLE', 'M13_SYSTSLCE', 'M13_JUVERHEU|M13_JUVEARTH'),
  function(x, y) {
    tmp1 <- filter(karnes, grepl(x, phewas_string, ignore.case=T), p<1e-3)[, c('snp', 'odds_ratio', 'p')]
    tmp2 <- filter(results.anno.stats.fdr.replicated, grepl(y, pheno, ignore.case=T))[, c('locus', 'beta', 'p.value')] 
    tmp2$beta <- exp(tmp2$beta)
    colnames(tmp2) <- colnames(tmp1)
    tmp <- inner_join(tmp1, tmp2, by='snp') %>% na.omit
    colnames(tmp) <- c('snp', 'odds_ratio_study', 'p_study', 'odds_ratio_y', 'p_y')
    if(nrow(tmp)>0) { data.frame(tmp, phenotype=x) } else NA
  })
karnes.comp <- karnes.comp[!is.na(karnes.comp)] %>% do.call(rbind, .) %>% distinct
karnes.comp$phenotype <- gsub('\\$|\\^', '', karnes.comp$phenotype)
karnes.comp$study <- 'Karnes et al. (2017)'

liu.comp <- map2(
  c('diabetes with other specified manifestations', 'celiac disease', 'rheumatoid arthritis', 'psoriasis$', 'ketoacidosis', 
    'pneumonia', 'systemic sclerosis', 'lichen planus'),
  c('T1D', 'K11_COELIAC', 'M13_RHEUMA', 'psoriasis', 'ketoacidosis',
    'pneumonia', 'M13_SYSTSLCE', 'L12_LICHENPLANUS'),
  function(x, y) {
    tmp1 <- filter(liu, grepl(x, Phenotype, ignore.case=T), p<1e-3)[, c('Gene', 'OR', 'p')]
    tmp2 <- filter(results.anno.stats.fdr.replicated, grepl(y, pheno, ignore.case=T))[, c('locus', 'beta', 'p.value')] 
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

hirata.comp <- map2(
  c('T1D', 'nephrotic syndrome', 'rheumatoid arthritis', 'asthma'),
  c('T1D', 'N14_NEPHROTICSYND', 'M13_RHEUMA', 'J10_ASTHMA'),
  function(x, y) {
    tmp1 <- filter(hirata, grepl(x, pheno, ignore.case=T), pval<1e-3)[, c('gene', 'Effect', 'pval')]
    tmp2 <- filter(results.anno.stats.fdr.replicated, grepl(y, pheno, ignore.case=T))[, c('locus', 'beta', 'p.value')] 
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


studies.comp <- rbind(karnes.comp, liu.comp, hirata.comp)
studies.comp$phenotype <- gsub('type 1 ', '', studies.comp$phenotype)
studies.comp <- unite(studies.comp, group, phenotype, study, sep='_', remove=F)

# plot
p.studies.comp <- ggplot(studies.comp, aes(odds_ratio_study, odds_ratio_y)) +
  geom_point(aes(fill=phenotype, shape=study, group=group, color=phenotype), alpha=0.95, size=4) +
  geom_label_repel(aes(odds_ratio_study, odds_ratio_y, label=snp), segment.size=0.15, nudge_x=0.5, nudge_y=-0.5, 
                   data=studies.comp[c(184, 187, 189, 54, 66, 109, 194, 199, 206, 190, 193, 164, 204, 208), ],
                   size=2.3, segment.alpha=0.6, color='grey25', alpha=0.8, force=8, min.segment.length=0, direction='both') +
  geom_label_repel(aes(odds_ratio_study, odds_ratio_y, label=snp), segment.size=0.15, nudge_x=-0.1, nudge_y=-0.2, 
                   data=studies.comp[c(52, 78, 138, 139, 140, 201), ],
                   size=2.3, segment.alpha=0.6, color='grey25', alpha=0.8, force=8, min.segment.length=0, direction='both') +
  geom_label_repel(aes(odds_ratio_study, odds_ratio_y, label=snp), segment.size=0.15, nudge_x=0.5, nudge_y=0,
                   data=studies.comp[c(116, 118, 85, 2, 125), ], 
                   direction='both', size=2.3, segment.alpha=0.6, color='grey25', alpha=0.8, force=8, min.segment.length=0) +
  scale_shape_manual(values=c(24, 21, 23)) +
  scale_color_manual(values=rep('white', studies.comp$phenotype %>% unique %>% length)) +
  scale_fill_manual(guide=guide_legend(override.aes=list(shape=22)),
    values=c(brewer.pal(8, name='Set2'), brewer.pal(8, name='Dark2'))[1:(studies.comp$phenotype %>% unique %>% length)]) +
  xlab('odds ratio (previous studies)') + ylab('odds ratio (present study)') +
  xlim(c(0, 7.1)) + ylim(c(0, 6)) +
  geom_hline(yintercept=1, size=0.2, linetype='dashed', color='grey40') +
  geom_vline(xintercept=1, size=0.2, linetype='dashed', color='grey40') +
  theme(panel.background=element_blank(),  legend.key=element_blank(),
        axis.text.x=element_text(size=8, angle=0, hjust=1), axis.text.y=element_text(size=8.5), 
        axis.line.x=element_line(size=.3), axis.line.y=element_line(size=.3), 
        axis.title.y=element_text(color='black'), legend.position='right',
        strip.background=element_rect(fill=NA),
        legend.key.size=unit(0.6, 'cm'), legend.text=element_text(size=8)) 

p.studies.comp.2 <- ggplot(studies.comp, aes(odds_ratio_study %>% log, odds_ratio_y %>% log)) +
  geom_point(aes(fill=phenotype, shape=study, group=group, color=phenotype), alpha=0.95, size=3) +
  geom_label_repel(aes(odds_ratio_study %>% log, odds_ratio_y %>% log, label=snp), segment.size=0.15, nudge_x=0.7, nudge_y=-0.5, 
                   data=studies.comp[c(184, 187, 189, 54, 66, 109, 194, 199, 206, 190, 193, 164, 204, 208), ],
                   size=1.8, segment.alpha=0.6, color='grey25', alpha=0.8, force=8, min.segment.length=0, direction='both') +
  geom_label_repel(aes(odds_ratio_study %>% log, odds_ratio_y %>% log, label=snp), segment.size=0.15, nudge_x=-0.2, nudge_y=0.5, 
                   data=studies.comp[c(52, 78, 138, 139, 140, 201), ],
                   size=1.8, segment.alpha=0.6, color='grey25', alpha=0.8, force=8, min.segment.length=0, direction='both') +
  geom_label_repel(aes(odds_ratio_study %>% log, odds_ratio_y %>% log, label=snp), segment.size=0.1, nudge_x=0.5, nudge_y=0,
                   data=studies.comp[c(116, 118, 85, 2, 125), ], 
                   direction='both', size=1.8, segment.alpha=0.6, color='grey25', alpha=0.8, force=8, min.segment.length=0) +
  scale_shape_manual(values=c(24, 21, 23)) +
  scale_color_manual(values=rep('white', studies.comp$phenotype %>% unique %>% length)) +
  scale_fill_manual(guide=guide_legend(override.aes=list(shape=22)),
                    values=c(brewer.pal(8, name='Set2'), brewer.pal(8, name='Dark2'))[1:(studies.comp$phenotype %>% unique %>% length)]) +
  xlab('beta (previous studies)') + ylab('beta (present study)') +
  xlim(c(-2, 2.5)) + ylim(c(-2, 2)) +
  geom_hline(yintercept=0, size=0.2, linetype='dashed', color='grey40') +
  geom_vline(xintercept=0, size=0.2, linetype='dashed', color='grey40') +
  theme(panel.background=element_blank(),  legend.key=element_blank(),
        axis.text.x=element_text(size=8, angle=0, hjust=1), axis.text.y=element_text(size=8.5), 
        axis.line.x=element_line(size=.3), axis.line.y=element_line(size=.3), 
        axis.title.y=element_text(color='black'), legend.position='right',
        strip.background=element_rect(fill=NA),
        legend.key.size=unit(0.6, 'cm'), legend.text=element_text(size=8)) 

## --------------------------------------------------------
## replication correlation
## --------------------------------------------------------

p.rep.corr <-ggplot(rbind(data.frame(filter(results.anno, beta<5.1, beta.replication<5.1)[c('beta', 'beta.replication')], group='all'),
                          data.frame(results.anno.stats.fdr.replicated[c('beta', 'beta.replication')], 
                                     group='FDR<0.01 & replicated')),  
                    aes(x=beta, y=beta.replication)) +
  geom_point(shape=21, color='grey10', alpha=0.4, size=0.3) +
  geom_smooth(method='lm', formula=y~x, size=0.3, alpha=0.7, se=F) +
  stat_cor(method="pearson", aes(label=..r.label..), size=3) +
  xlab('beta (discovery cohort)') + ylab('beta (replication cohort)') +
  facet_wrap(~group, scales='free', nrow=2, ncol=1) +
  theme(strip.background=element_rect(fill='grey90', color='white'), strip.text=element_text(size=10),
        panel.border=element_blank(), panel.background=element_blank(), axis.line.x=element_line(size=0.3),
        axis.line.y=element_line(size=0.3), axis.text.y=element_text(size=8.5), axis.text.x=element_text(size=8.5),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())



## --------------------------------------------------------
## plot of shared alleles in infections and autoimmune
## --------------------------------------------------------

# conditional analysis results
results.cond.f <- filter(results.cond) %>% na.omit#, p.value<5e-8)
results.cond.f$HLA <- formHLA(results.cond.f$HLA)
# write
write.table(results.cond.f, './results/tables/phenotype_conditional_23Jun20.tsv', quote=F, sep='\t', row.names=F)


# select infect/autoimm phenos
results.cond.f.infect <- lapply(results.cond.f$HLA %>% unique, function(i) {
  out <- filter(results.cond.f, HLA==i)
  out <- out[grepl('AB1|INFECT|ITIS|SEPSIS', out$pheno), ]
  out <- out[grepl('G6_DEMYEL|E4_DIABETES|J10_ASTHMA|L12_PSORIASIS|M13_ANKYLOSPON|M13_RHEUMA', out$covpheno), ]
  return(out)
}) %>% do.call(rbind, .)
results.cond.f.autoimm <- lapply(results.cond.f$HLA %>% unique, function(i) {
  out <- filter(results.cond.f, HLA==i)
  out <- out[grepl('G6_DEMYEL|E4_DIABETES|J10_ASTHMA|L12_PSORIASIS|M13_ANKYLOSPON|M13_RHEUMA', out$pheno), ]
  out <- out[grepl('AB1|INFECT|ITIS|SEPSIS', out$covpheno), ]
  return(out)
}) %>% do.call(rbind, .)

# select HLAs that independently associate to both infect. and autoimm
results.cond.f.infect.2 <- filter(results.cond.f.infect, HLA %in% intersect(results.cond.f.autoimm$HLA, results.cond.f.infect$HLA))
results.cond.f.autoimm.2 <- filter(results.cond.f.autoimm, HLA %in% intersect(results.cond.f.autoimm$HLA, results.cond.f.infect$HLA))

# plot
p.infect.cond <- ggplot(results.cond.f.infect.2, aes(HLA, covpheno, fill=beta)) +
  geom_tile(color='black') +
  ylab('covariate phenotype') +
  facet_wrap(~pheno) +
  scale_fill_gradient2(low="seagreen", high="red", na.value="grey") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7), axis.text.y=element_text(size=7.5), 
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0), size=12),
        axis.title.x=element_blank(),
        panel.border=element_rect(size=0.3, fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white", colour="white"), axis.ticks=element_blank(),
        legend.key.size=unit(0.75, 'line'), legend.margin=margin(t=0, r=10, b=0, l=10),
        legend.title=element_text(size=10), legend.text=element_text(size=8))  

p.autoimm.cond <- ggplot(results.cond.f.autoimm.2, aes(HLA, covpheno, fill=beta)) +
  geom_tile(color='black') +
  ylab('covariate phenotype') +
  facet_wrap(~pheno) +
  scale_fill_gradient2(low="seagreen", high="red", na.value="grey") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7), axis.text.y=element_text(size=7.5), 
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0), size=12),
        axis.title.x=element_blank(),
        panel.border=element_rect(size=0.3, fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.key.size=unit(0.75, 'line'), legend.margin=margin(t=0, r=10, b=0, l=10),
        panel.background=element_rect(fill="white", colour="white"), axis.ticks=element_blank(), 
        legend.title=element_text(size=10), legend.text=element_text(size=8)) 


## --------------------------------------------------------
## pleiotropy ccFDR 
## --------------------------------------------------------

# https://rdrr.io/github/alexploner/condFDR/man/ccFDR.html

tmp.infectious <- results.cond.f.infect.2$pheno %>% unique
tmp.autoimmune <- results.cond.f.autoimm.2$pheno %>% unique

inf.autoimm.ccfdr <- map(tmp.infectious, function(i) {
  map(tmp.autoimmune, function(a) {
    rep  <- filter(results.anno.stats.fdr.replicated, pheno==i | pheno==a)$locus
    out1 <- filter(results.anno, pheno==i, locus %in% rep)
    out2 <- filter(results.anno, pheno==a, locus %in% rep)
    out  <- inner_join(out1, out2, by='locus')
    out  <- ccFDR(out, p1='FisherP.x', p2='FisherP.y')[, c('locus', 'ccFDR')]
    out$phenoi <- i
    out$phenoa <- a
    return(out)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)
expand.grid(tmp.infectious, tmp.autoimmune) %>% nrow
inf.autoimm.ccfdr$adj_ccFDR <- inf.autoimm.ccfdr$ccFDR*35
inf.autoimm.ccfdr$log_ccFDR <- -log10(inf.autoimm.ccfdr$adj_ccFDR)

p.inf.autoimm.ccfdr <- ggplot(filter(inf.autoimm.ccfdr, adj_ccFDR<0.001), aes(locus, phenoa, fill=log_ccFDR)) +
  geom_tile(color='black') +
  facet_wrap(~phenoi) +
  scale_fill_gradient2(low="white", high="red", na.value="grey", name=expression('-log'[10]*'(ccFDR)')) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7), axis.text.y=element_text(size=7.5), 
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0), size=12),
        axis.title.x=element_blank(),
        panel.border=element_rect(size=0.3, fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.key.size=unit(0.75, 'line'), legend.margin=margin(t=0, r=10, b=0, l=10),
        panel.background=element_rect(fill="white", colour="white"), axis.ticks=element_blank(), 
        legend.title=element_text(size=10), legend.text=element_text(size=8), legend.position='top') 

## --------------------------------------------------------
## combine ccFDR and infections and autoimmune
## conditional analysis
## --------------------------------------------------------

inf.autoimm.ccfdr
results.cond.f.infect.2
results.cond.f.autoimm.2

results.cond.f.infect.3 <- inner_join(results.cond.f.infect.2, inf.autoimm.ccfdr, 
                                      by=c('pheno'='phenoi', 'covpheno'='phenoa', 'HLA'='locus'))
results.cond.f.autoimm.3 <- inner_join(results.cond.f.autoimm.2, inf.autoimm.ccfdr, 
                                      by=c('pheno'='phenoa', 'covpheno'='phenoi', 'HLA'='locus'))

# plot
p.infect.cond <- ggplot(results.cond.f.infect.3 %>% filter(ccFDR<0.0002857143), aes(HLA, covpheno, fill=beta)) +
  geom_tile(color='black') +
  ylab('covariate phenotype') +
  facet_wrap(~pheno) +
  scale_fill_gradient2(low="seagreen", high="red", na.value="grey") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7), axis.text.y=element_text(size=7.5), 
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0), size=12),
        axis.title.x=element_blank(),
        panel.border=element_rect(size=0.3, fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white", colour="white"), axis.ticks=element_blank(),
        legend.key.size=unit(0.75, 'line'), legend.margin=margin(t=0, r=10, b=0, l=10),
        legend.title=element_text(size=10), legend.text=element_text(size=8))  

p.autoimm.cond <- ggplot(results.cond.f.autoimm.3%>% filter(ccFDR<0.0002857143), aes(HLA, covpheno, fill=beta)) +
  geom_tile(color='black') +
  ylab('covariate phenotype') +
  facet_wrap(~pheno) +
  scale_fill_gradient2(low="seagreen", high="red", na.value="grey") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7), axis.text.y=element_text(size=7.5), 
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0), size=12),
        axis.title.x=element_blank(),
        panel.border=element_rect(size=0.3, fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.key.size=unit(0.75, 'line'), legend.margin=margin(t=0, r=10, b=0, l=10),
        panel.background=element_rect(fill="white", colour="white"), axis.ticks=element_blank(), 
        legend.title=element_text(size=10), legend.text=element_text(size=8)) 




## --------------------------------------------------------
## Figures
## --------------------------------------------------------

# Figure S1

# Figure 1
jpeg('./results/figures/Figure1.jpg', height=4.5, width=10, units='in', res=600)                           
p.hla.phewas
dev.off()

# Figure 2
jpeg('./results/figures/Figure2.jpg', height=5.5, width=10, units='in', res=600)                           
ggarrange(p.studies.comp.2, p.rep.corr, 
          nrow=1, ncol=2, widths=c(2, 0.85), labels=c('a', 'b'), font.label=list(size=15))
dev.off()

# Figure 3
jpeg('./results/figures/Figure3.jpg', height=7, width=10, units='in', res=600)                           
ggarrange(p.infect.cond, p.autoimm.cond, 
          nrow=2, ncol=1, heights=c(1.05, 0.95), labels=c('a', 'b'), common.legend=T, align='v', 
          legend='right', font.label=list(size=15))
dev.off()


