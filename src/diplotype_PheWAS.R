

## Analysis of HLA diplotype PheWAS results

source('./src/functions.R')


## --------------------------------------------------------
## data
## --------------------------------------------------------

# read resukts and clean
geno.results <- fread('./data/R3_genotype_phewas_cases_over4_30Mar20.tsv', data.table=F)
geno.results <- geno.results[str_split_fixed(geno.results$risk.allele, ':', 2)[, 2]!='', ] # keep only 4-digit
geno.results <- geno.results[str_split_fixed(geno.results$mod.allele, ':', 2)[, 2]!='', ] # keep only 4-digit
geno.results$risk.allele <- gsub('DQB|DRB|DQA|DPB', '', geno.results$risk.allele) # clean
geno.results$mod.allele  <- gsub('DQB|DRB|DQA|DPB', '', geno.results$mod.allele) # clean
geno.results <- filter(geno.results, !is.na(p.value))

# replication results
geno.results.rep <- fread('./data/R5sansR3_genotype_phewas_cases_over4_25Mar20.tsv', data.table=F)
geno.results.rep <- geno.results.rep[str_split_fixed(geno.results.rep$risk.allele, ':', 2)[, 2]!='', ] # keep only 4-digit
geno.results.rep <- geno.results.rep[str_split_fixed(geno.results.rep$mod.allele, ':', 2)[, 2]!='', ] # keep only 4-digit
geno.results.rep$risk.allele <- gsub('DQB|DRB|DQA|DPB', '', geno.results.rep$risk.allele) # clean
geno.results.rep$mod.allele  <- gsub('DQB|DRB|DQA|DPB', '', geno.results.rep$mod.allele) # clean
geno.results.rep <- filter(geno.results.rep, !is.na(p.value))

# fdr
geno.results     <- mutate(geno.results, adjusted.p=adaptiveBH(p.value, alpha=0.01, silent=T)$adjPValues)
geno.results.fdr <- filter(geno.results, adjusted.p<0.01)

# stats
geno.results.stats     <- fread('./data/R3_geno_stats_n5.tsv', data.table=F)
geno.results.rep.stats <- fread('./data/R5_geno_stats_n5.tsv', data.table=F)


## --------------------------------------------------------
## phenotype annotation 
## --------------------------------------------------------

# read phenotype annotation
fg.anno      <- fread('~/Documents/FinnGen/FINNGEN_ENDPOINTS_DF5_V2_2020-02-11_public.tsv', data.table=F)[, c(1,4,5)]
fg.anno$TAGS <- gsub('#', '', fg.anno$TAGS, fixed=T)
fg.anno$TAGS <- str_split_fixed(fg.anno$TAGS, ',', 10)[, 1]

# join annotation to data
geno.results.anno <- geno.results[, c(1:5, 14, 8:13)]
geno.results.anno <- left_join(geno.results.anno, fg.anno, by=c('pheno'='NAME'))
geno.results.anno <- geno.results.anno[, c(13, 1, 14, 2:12)]
colnames(geno.results.anno)[c(1:3)] <- c('tag', 'pheno', 'longname') 
geno.results.anno[is.na(geno.results.anno)] <- ''

# join imputation stats to data
geno.results.anno  <- unite(geno.results.anno, phenoallele, pheno, locus, risk.allele, mod.allele, sep='_', remove=F)
geno.results.stats <- mutate(geno.results.stats, risk.allele=str_split_fixed(geno.results.stats$locus, '\\.', 2)[, 2])
geno.results.stats$locus      <- str_split_fixed(geno.results.stats$locus, '\\.', 2)[, 1]
geno.results.stats$mod.allele <- gsub('.', ':', geno.results.stats$mod.allele, fixed=T)
geno.results.stats <- unite(geno.results.stats, phenoallele, pheno, locus, risk.allele, mod.allele, sep='_', remove=F)
geno.results.anno  <- left_join(geno.results.anno, geno.results.stats[, c(1, 5:6)], by='phenoallele')

# include replication results
geno.results.rep2d <- unite(geno.results.rep, phenoallele, pheno, locus, risk.allele, mod.allele, sep='_') %>% 
  .[, c('phenoallele', 'p.value', 'beta', 'SEbeta')]
colnames(geno.results.rep2d)[2:4] <- c('p.value.replication', 'beta.replication', 'SEbeta.replication')
geno.results.anno <- left_join(geno.results.anno, geno.results.rep2d, by='phenoallele')

# include replication imputation stats
geno.results.rep.stats <- mutate(geno.results.rep.stats, risk.allele=str_split_fixed(geno.results.rep.stats$locus, '\\.', 2)[, 2])
geno.results.rep.stats$locus      <- str_split_fixed(geno.results.rep.stats$locus, '\\.', 2)[, 1]
geno.results.rep.stats$mod.allele <- gsub('.', ':', geno.results.rep.stats$mod.allele, fixed=T)
geno.results.rep.stats$mod.allele <- str_pad(geno.results.rep.stats$mod.allele, 5, 'left', pad='0')
geno.results.rep.stats <- unite(geno.results.rep.stats, phenoallele, pheno, locus, risk.allele, mod.allele, sep='_', remove=F)
colnames(geno.results.rep.stats)[5:12] <- paste0(colnames(geno.results.rep.stats)[5:12], '.replication')
geno.results.anno <- left_join(geno.results.anno, geno.results.rep.stats[, c(1, 7:8, 10:11, 5:6)], by='phenoallele')
colnames(geno.results.anno) <- gsub('geno.cases.replication', 'geno.n.cases.replication', colnames(geno.results.anno), fixed=T)

# remove duplicated if exists
geno.results.d    <- unite(geno.results.anno, phenoallele, pheno, locus, risk.allele, mod.allele, sep='_')
geno.results.anno <- geno.results.anno[!duplicated(geno.results.d$phenoallele), ] %>% dplyr::select(-phenoallele)

# replicated
geno.results.anno.replicated <- filter(geno.results.anno, adjusted.p<0.01 & p.value.replication<0.01)

# write out
geno.results.anno.replicated.list <- tapply(1:nrow(geno.results.anno.replicated), geno.results.anno.replicated$pheno, 
                                            function(x) geno.results.anno.replicated[x, ]) 
write.table(lapply(geno.results.anno.replicated.list, function(x) 
  rbind(data.frame(x), rep(' ', ncol(geno.results.anno.replicated)))) %>% 
              do.call(rbind, .), './results/tables/PheWAS_genotype_replicated_12June20.tsv', sep='\t', row.names=F)


## --------------------------------------------------------
## test beta differences
## --------------------------------------------------------

# compute z test for beta differences
geno.results.f <- na.omit(geno.results.anno)
geno.results.z <- tapply(1:nrow(geno.results.f), geno.results.f$pheno, function(x) {
  tmp <- geno.results.f[x, ] %>% arrange(desc(beta))
  #tmp <- filter(tmp, imp.pprob.mean>0.6, imp.pprob.mean.replication>0.6)
  out <- sapply(2:nrow(tmp), function(i) {
    Z.test.2(tmp[1, 'beta'], tmp[i, 'beta'], tmp[1, 'SEbeta'], tmp[i, 'SEbeta'])
  })
  data.frame(tmp, beta.diff.p.value=c(NA, out))
}) %>% do.call(rbind, .)
geno.results.z <- mutate(geno.results.z, adjusted.beta.diff.p.value=adaptiveBH(beta.diff.p.value, alpha=0.01, silent=T)$adjPValues)
geno.results.z <- geno.results.z[
  ifelse(geno.results.z$adjusted.beta.diff.p.value<0.01 | is.na(geno.results.z$adjusted.beta.diff.p.value), T, F), ]
geno.results.z <- geno.results.z[!duplicated(geno.results.z), ]
geno.results.z <- geno.results.z %>% add_count(pheno) %>% filter(n>1) %>% dplyr::select(-n) %>% data.frame

# keep phenos that have a significant beta diff
geno.results.z <- lapply(geno.results.z$pheno %>% unique, function(x) {
  tmp <- filter(geno.results.z, pheno==x)
  out <- NA
  if(any(na.omit(tmp$adjusted.beta.diff.p.value)<0.01)) out <- tmp else out <- NA
  out
}) 
geno.results.z <- geno.results.z[!is.na(geno.results.z)] %>% do.call(rbind, .)
geno.results.z <- geno.results.z[!duplicated(geno.results.z), ]
geno.results.z <- geno.results.z %>% add_count(pheno) %>% filter(n>1) %>% dplyr::select(-n) %>% data.frame

# beta diff test for replication
geno.results.z <- tapply(1:nrow(geno.results.z), geno.results.z$pheno, function(x) {
  tmp <- geno.results.z[x, ] #%>% arrange(desc(beta))
  out <- sapply(2:nrow(tmp), function(i) {
    Z.test.2(tmp[1, 'beta.replication'], tmp[i, 'beta.replication'], tmp[1, 'SEbeta.replication'], tmp[i, 'SEbeta.replication'])
  })
  data.frame(tmp, replication.beta.diff.p.value=c(NA, out))
}) %>% do.call(rbind, .)

keep.phe <- tapply(1:nrow(geno.results.z), geno.results.z$pheno, function(x) {
  tmp <- geno.results.z[x, ]
  tmp$adjusted.p[1]<0.01 & tmp$p.value.replication[1]<0.01
}) 
keep.phe <- keep.phe[keep.phe==T] %>% names
geno.results.z <- geno.results.z[geno.results.z$pheno %in% keep.phe, ]

# check beta diffs are replicated
geno.results.z <- filter(geno.results.z, replication.beta.diff.p.value<0.01 | is.na(replication.beta.diff.p.value))
geno.results.z <- geno.results.z %>% add_count(pheno) %>% filter(n>1) %>% dplyr::select(-n) %>% data.frame


# filter by difference metric between discovery and replication
geno.results.z.clean <- mutate(geno.results.z, betarepdiff=(beta-beta.replication))
geno.results.z.clean <- filter(geno.results.z.clean, 
                               betarepdiff < (mean(betarepdiff) + 1.5*sd(betarepdiff)), 
                               betarepdiff > (mean(betarepdiff) - 1.5*sd(betarepdiff)))

# write out
write.table(geno.results.z, './results/R3_genotype_betadiff_replicated_18May20.tsv', sep='\t', row.names=F)

geno.results.z.list <- tapply(1:nrow(geno.results.z), 
                              geno.results.z$pheno, function(x) geno.results.z[x, ]) 
write.table(lapply(geno.results.z.list, function(x) 
  rbind(data.frame(x), rep(' ', ncol(geno.results.z)))) %>% 
    do.call(rbind, .), './results/tables/PheWAS_genotype_betadiffs_replicated_15Jun20.tsv', sep='\t', row.names=F)


## --------------------------------------------------------
## beta distribution 
## --------------------------------------------------------

geno.results.z.diffs <- (geno.results.z$pheno %>% table %>% names) %>% map(function(x) {
  tt       <- filter(geno.results.z, pheno==x)#, imp.pprob.mean>0.6, imp.pprob.mean.replication>0.6)
  ind.min  <- which.min(tt$beta)
  ind.max  <- which.max(tt$beta)
  beta.min <- tt[ind.min, c('beta', 'beta.replication')] %>% rowMeans
  beta.max <- tt[ind.max, c('beta', 'beta.replication')] %>% rowMeans
  p.rep    <- tt[ind.min, 'replication.beta.diff.p.value'] 
  data.frame(pheno=x, betadiff=abs(beta.min-beta.max), p.rep=p.rep)
}) %>% do.call(rbind, .)
geno.results.z.diffs <- arrange(geno.results.z.diffs, desc(betadiff))
geno.results.z.diffs <- filter(geno.results.z.diffs, p.rep<0.01)


# summary figure of interactions vs allele-level
geno.results.z$locus %>% table
filter(geno.results.z, locus=='DQB1')$risk.allele %>% table
filter(geno.results.z, locus=='DQB1')$mod.allele %>% table
filter(geno.results.z, locus=='DQB1', risk.allele=='03:02')$tag %>% table
filter(geno.results.z, locus=='DRB1')$risk.allele %>% table
filter(geno.results.z, locus=='DQA1')$risk.allele %>% table
filter(geno.results.z, locus=='B')$risk.allele %>% table

alleleVsDiplo <- function(gene, allele) {
  
  genallel <- paste0(gene, '*', allele)
  out <- filter(geno.results.z, locus==gene, risk.allele==allele)
  ttt <- filter(results.anno.stats, pheno %in% unique(out$pheno))
  out <- filter(out, pheno %in% ttt$pheno)
  
  out <- map(out$pheno %>% unique, function(x) {
    
    tmp01 <- filter(ttt, pheno==x)
    genallel.beta <- filter(tmp01, locus==genallel)[, c('beta', 'beta.replication')] %>% rowMeans()
    genallel.sebeta <- filter(tmp01, locus==genallel)[, c('SEbeta', 'SEbeta.replication')] %>% rowMeans()
    tmp02 <- filter(out, pheno==x)
    tmp01 <- tmp01[match(paste(tmp02$locus, tmp02$mod.allele, sep='*'), tmp01$locus), ]
    
    tmp01 <- data.frame(tmp01[, c('tag', 'pheno', 'longname', 'locus')], 
                        risk.gene=genallel,
                        risk.beta=genallel.beta, 
                        risk.sebeta=genallel.sebeta, 
                        mbeta=rowMeans(tmp01[, c('beta', 'beta.replication')], na.rm=T),
                        sembeta=rowMeans(tmp01[, c('SEbeta', 'SEbeta.replication')], na.rm=T))
    tmp01$comp.beta <- tmp01$risk.beta + tmp01$mbeta
    tmp01$comp.sebeta <- tmp01$risk.sebeta + tmp01$sembeta
    data.frame(tmp01, 
               ipheno=tmp02[, c('pheno')], 
               risk.allele=tmp02$risk.allele,
               mod.allele=tmp02$mod.allele,
               interact.beta=rowMeans(tmp02[, c('beta', 'beta.replication')]),
               interact.sebeta=rowMeans(tmp02[, c('SEbeta', 'SEbeta.replication')]))
    
  }) %>% do.call(rbind, .)
  out[!is.na(out$mbeta), ]
  
}


geno.betacomparison <- map2(c('DQB1',  'DQB1',  'DQB1',  'DRB1',  'DRB1',  'DRB1',  'B',     'DQA1'), 
                            c('03:02', '02:01', '05:01', '04:01', '04:08', '15:01', '27:05', '05:01'), 
                            function(x, y) { alleleVsDiplo(x ,y) }) %>% do.call(rbind, .)
geno.betacomparison$tag[grepl('DM|DIAB|Diab|INSULIN|E4', geno.betacomparison$pheno)] <- 'Diabetes'
geno.betacomparison$tag[grepl('DM|DIAB|Diab', geno.betacomparison$tag)] <- 'Diabetes'
geno.betacomparison$tag[geno.betacomparison$pheno=='AUTOIMMUNE'] <- 'Autoimmune'
geno.betacomparison$tag[geno.betacomparison$tag=='M13'] <- 'Rheumatoid arthritis'
geno.betacomparison$tag[geno.betacomparison$tag=='K11'] <- 'Diseases of the digestive system'
geno.betacomparison$tag[geno.betacomparison$tag=='H7'] <- 'Diseases of the eye'
geno.betacomparison$tag[geno.betacomparison$tag=='ILD_CM'] <- 'ILD-related Thyroiditis'
geno.betacomparison$tag[geno.betacomparison$tag=='RHEUMA'] <- 'Rheumatoid arthritis'
geno.betacomparison$tag[geno.betacomparison$tag=='RHEUMA_CM'] <- 'Rheumatoid arthritis'
geno.betacomparison$tag[geno.betacomparison$tag=='N14'] <- 'Glomerular diseases'
geno.betacomparison <- filter(geno.betacomparison, tag!='RX') %>% filter(., tag!='Z21')
geno.betacomparison$tag2 <- factor(geno.betacomparison$tag, 
                                  levels=c(geno.betacomparison$tag, haplo.betacomparison$tag) %>% unique)
geno.betacomparison$diff <- geno.betacomparison$interact.beta - geno.betacomparison$comp.beta


geno.betacomparison.f <- map(geno.betacomparison$pheno %>% unique, function(x) {
  out <- filter(geno.betacomparison, pheno==x)
  arrange(out, 1/abs(diff))[1, c(1:2, 4,3, 5:7, 10:11, 15:16, 18)]
}) %>% do.call(rbind, .) %>% arrange(., 1/(abs(diff))) %>% 
  filter(., pheno %in% c('K11_COELIAC', 'T1D_STRICT1', 'RHEUMA_SEROPOS', 'G6_DEMYEL' , 'JUVEN_ARTHR_COMORB', 
                         'H7_IRIDOCYCLITIS', 'N14_GLOMEINOTH', 'THYROIDITIS_ILD')) %>% arrange(., diff)
geno.betacomparison.f$pheno <- paste0('(', geno.betacomparison.f$pheno ,')')
geno.betacomparison.f <- unite(geno.betacomparison.f, fullname, longname,  pheno, sep=' ')
geno.betacomparison.f <- unite(geno.betacomparison.f, HLAs, risk.gene, locus, sep=' - ')
geno.betacomparison.f$fullname <- factor(geno.betacomparison.f$fullname, levels=geno.betacomparison.f$fullname %>% unique)

