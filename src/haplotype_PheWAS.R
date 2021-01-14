

## Analysis of HLA haplotype PheWAS results

source('./src/functions.R')


## --------------------------------------------------------
## data 
## --------------------------------------------------------

# read and clean
haplo.results <- fread('./data/R3_haplo_phewas_casesn_over4_25Mar20.tsv', data.table=F)
haplo.results <- haplo.results[str_split_fixed(haplo.results$risk.locus, ':', 2)[, 2]!='', ] # keep only 4-digit
haplo.results <- haplo.results[str_split_fixed(haplo.results$mod.locus, '.', 3)[, 3]!='', ] # keep only 4-digit
haplo.results <- filter(haplo.results, !is.na(p.value))
colnames(haplo.results)[1] <- 'pheno'
haplo.results$risk.locus <- formHLA(haplo.results$risk.locus)
haplo.results$mod.locus <- formHLA(haplo.results$mod.locus)

# replication
haplo.results.rep <- fread('./data/R5sansR3_haplo_phewas_casesn_over4_25Mar20.tsv', data.table=F)
haplo.results.rep <- haplo.results.rep[str_split_fixed(haplo.results.rep$risk.locus, ':', 2)[, 2]!='', ] # keep only 4-digit
haplo.results.rep <- haplo.results.rep[str_split_fixed(haplo.results.rep$mod.locus, '.', 3)[, 3]!='', ] # keep only 4-digit
haplo.results.rep <- filter(haplo.results.rep, !is.na(p.value))
colnames(haplo.results.rep)[1] <- 'pheno'
haplo.results.rep$risk.locus <- formHLA(haplo.results.rep$risk.locus)
haplo.results.rep$mod.locus <- formHLA(haplo.results.rep$mod.locus)

# fdr
haplo.results     <- mutate(haplo.results, adjusted.p=adaptiveBH(p.value, alpha=0.01, silent=T)$adjPValues) 
haplo.results.fdr <- haplo.results %>% filter(., adjusted.p<0.01)

# stats
haplo.results.stats     <- fread('./data/R3_haplo_stats_n5.tsv', data.table=F)
haplo.results.rep.stats <- fread('./data/R5_haplo_stats_n5.tsv', data.table=F)


## --------------------------------------------------------
## annotate phenos
## --------------------------------------------------------

# phenotype annotation
fg.anno      <- fread('~/Documents/FinnGen/FINNGEN_ENDPOINTS_DF5_V2_2020-02-11_public.tsv', data.table=F)[, c(1,4,5)]
fg.anno$TAGS <- gsub('#', '', fg.anno$TAGS, fixed=T)
fg.anno$TAGS <- str_split_fixed(fg.anno$TAGS, ',', 10)[, 1]

# join annotation to results
haplo.results.anno <- left_join(haplo.results, fg.anno, by=c('pheno'='NAME'))
haplo.results.anno <- haplo.results.anno[, c(14, 1, 15, 2:4, 13, 7:12)]
colnames(haplo.results.anno)[c(1:3)] <- c('tag', 'pheno', 'longname') 

# join to imputation stats
colnames(haplo.results.stats) <- gsub('allele', 'locus', colnames(haplo.results.stats))
colnames(haplo.results.stats) <- gsub('haplo.', '', colnames(haplo.results.stats), fixed=T)
haplo.results.stats$mod.locus <- formHLA(haplo.results.stats$mod.locus)
haplo.results.stats           <- unite(haplo.results.stats, phenoallele, pheno, risk.locus, mod.locus, sep='_', remove=F)
haplo.results.anno  <- unite(haplo.results.anno, phenoallele, pheno, risk.locus, mod.locus, sep='_', remove=F)
haplo.results.anno  <- left_join(haplo.results.anno, haplo.results.stats[, c(1, 7:10)], by='phenoallele')

# include replication results
colnames(haplo.results.rep) <- paste0(colnames(haplo.results.rep), '.replication')
haplo.results.rep           <- unite(haplo.results.rep, phenoallele, pheno.replication, 
                                     risk.locus.replication, mod.locus.replication, sep='_', remove=F)
haplo.results.anno          <- left_join(haplo.results.anno, haplo.results.rep[, c(1, 5, 8:13)], by='phenoallele')
  
# join with replication stats
colnames(haplo.results.rep.stats) <- gsub('allele', 'locus', colnames(haplo.results.rep.stats))
colnames(haplo.results.rep.stats) <- gsub('haplo.', '', colnames(haplo.results.rep.stats), fixed=T)
haplo.results.rep.stats$mod.locus <- formHLA(haplo.results.rep.stats$mod.locus)
colnames(haplo.results.rep.stats) <- paste0(colnames(haplo.results.rep.stats), '.replication')
haplo.results.rep.stats <- unite(haplo.results.rep.stats, phenoallele, pheno.replication, 
                                 risk.locus.replication, mod.locus.replication, sep='_', remove=F)
haplo.results.anno <- left_join(haplo.results.anno, haplo.results.rep.stats[, c(1, 7:10)], by='phenoallele')
colnames(haplo.results.anno) <- gsub('pp', 'imp.pprob', colnames(haplo.results.anno))
haplo.results.anno <- dplyr::select(haplo.results.anno, -phenoallele)
haplo.results.anno.replicated <- filter(haplo.results.anno, adjusted.p<0.01 & p.value.replication<0.01)

# write out
haplo.results.anno.list <- tapply(1:nrow(haplo.results.anno), haplo.results.anno$pheno, function(x) haplo.results.anno[x, ]) 
write.table(lapply(haplo.results.anno.list, function(x) rbind(data.frame(x), rep(' ', ncol(haplo.results.anno)))) %>% 
              do.call(rbind, .), './results/R3_haplotype_PheWAS_full_09June20.tsv', sep='\t', row.names=F)
haplo.results.anno.replicated.list <- tapply(1:nrow(haplo.results.anno.replicated), 
                                             haplo.results.anno.replicated$pheno, function(x) haplo.results.anno.replicated[x, ]) 
write.table(lapply(haplo.results.anno.replicated.list, function(x) 
  rbind(data.frame(x), rep(' ', ncol(haplo.results.anno.replicated)))) %>% 
              do.call(rbind, .), './results/tables/PheWAS_haplotype_replicated_09June20.tsv', sep='\t', row.names=F)


## --------------------------------------------------------
## calculate beta differences
## --------------------------------------------------------

haplo.results.f <- na.omit(haplo.results.anno)

# beta z-test within phenotypes
haplo.results.z <- tapply(1:nrow(haplo.results.f), haplo.results.f$pheno, function(x) {
  tmp <- haplo.results.f[x, ] %>% arrange(desc(beta))
  out <- sapply(2:nrow(tmp), function(i) {
    Z.test.2(tmp[1, 'beta'], tmp[i, 'beta'], tmp[1, 'SEbeta'], tmp[i, 'SEbeta'])
  })
  data.frame(tmp, beta.diff.p.value=c(NA, out))
}) %>% do.call(rbind, .)

# clean
haplo.results.z <- mutate(haplo.results.z, adjusted.beta.diff.p.value=adaptiveBH(beta.diff.p.value, alpha=0.01, silent=T)$adjPValues)
haplo.results.z <- haplo.results.z[
  ifelse(haplo.results.z$adjusted.beta.diff.p.value<0.01 | is.na(haplo.results.z$adjusted.beta.diff.p.value), T, F), ]
haplo.results.z <- haplo.results.z[!duplicated(haplo.results.z), ]
haplo.results.z <- haplo.results.z %>% add_count(pheno) %>% filter(n>1) %>% dplyr::select(-n) %>% data.frame

# beta diff test for replication
haplo.results.z <- tapply(1:nrow(haplo.results.z), haplo.results.z$pheno, function(x) {
  tmp <- haplo.results.z[x, ] %>% arrange(desc(beta))
  out <- sapply(2:nrow(tmp), function(i) {
    Z.test.2(tmp[1, 'beta.replication'], tmp[i, 'beta.replication'], tmp[1, 'SEbeta.replication'], tmp[i, 'SEbeta.replication'])
  })
  data.frame(tmp, replication.beta.diff.p.value=c(NA, out))
}) %>% do.call(rbind, .)

# check that best allele is associated and replicated
keep.phe <- tapply(1:nrow(haplo.results.z), haplo.results.z$pheno, function(x) {
  tmp <- haplo.results.z[x, ]
  tmp$adjusted.p[1]<0.01 & tmp$p.value.replication[1]<0.01
}) 
keep.phe <- keep.phe[keep.phe==T] %>% names
haplo.results.z <- haplo.results.z[haplo.results.z$pheno %in% keep.phe, ]

# check beta diffs are replicated
haplo.results.z <- filter(haplo.results.z, replication.beta.diff.p.value<0.01 | is.na(replication.beta.diff.p.value))
haplo.results.z <- haplo.results.z %>% add_count(pheno) %>% filter(n>1) %>% dplyr::select(-n) %>% data.frame

# clean for plotting
haplo.results.z.clean <- mutate(haplo.results.z, betarepdiff=(beta-beta.replication))
haplo.results.z.clean <- filter(haplo.results.z.clean, 
                               betarepdiff < (mean(betarepdiff) + 1.5*sd(betarepdiff)), 
                               betarepdiff > (mean(betarepdiff) - 1.5*sd(betarepdiff)))

# write out 
haplo.results.z.list <- tapply(1:nrow(haplo.results.z), haplo.results.z$pheno, function(x) haplo.results.z[x, ]) 
write.table(lapply(haplo.results.z.list, function(x) 
  rbind(data.frame(x), rep(' ', ncol(haplo.results.z.anno.replicated)))) %>% 
              do.call(rbind, .), './results/tables/PheWAS_haplotype_betadiffs_replication_09June20.tsv', sep='\t', row.names=F)


       
## --------------------------------------------------------
## beta distribution plot
## --------------------------------------------------------

# summary figure of interactions vs allele-level
haplo.results.z$risk.locus %>% table
filter(haplo.results.z, risk.locus==ć)$mod.locus %>% table
filter(haplo.results.z, risk.locus=='B*27:05')$mod.locus %>% table
filter(haplo.results.z, risk.locus=='DQB1*03:02')$mod.locus %>% table

alleleVshaplo <- function(gene1) {
  
  out <- filter(haplo.results.z, risk.locus==gene1)
  ttt <- filter(results.anno.stats, pheno %in% unique(out$pheno))
  out <- filter(out, pheno %in% ttt$pheno)
  
  out <- map(out$pheno %>% unique, function(x) {
    
    tmp01 <- filter(ttt, pheno==x)
    genallel.beta   <- filter(tmp01, locus==gene1)[, c('beta', 'beta.replication')] %>% rowMeans()
    genallel.sebeta <- filter(tmp01, locus==gene1)[, c('SEbeta', 'SEbeta.replication')] %>% rowMeans()
    tmp02 <- filter(out, pheno==x)
    tmp01 <- tmp01[match(tmp02$mod.locus, tmp01$locus), ]
    
    tmp01 <- data.frame(tmp01[, c('tag', 'pheno', 'longname', 'locus')], 
                        risk.gene=gene1,
                        risk.beta=genallel.beta, 
                        risk.sebeta=genallel.sebeta, 
                        mbeta=rowMeans(tmp01[, c('beta', 'beta.replication')], na.rm=T),
                        sembeta=rowMeans(tmp01[, c('SEbeta', 'SEbeta.replication')], na.rm=T))
    tmp01$comp.beta   <- tmp01$risk.beta + tmp01$mbeta
    tmp01$comp.sebeta <- tmp01$risk.sebeta + tmp01$sembeta
    data.frame(tmp01, 
               ipheno=tmp02[, c('pheno')], 
               risk.allele=tmp02$risk.locus %>% unique,
               mod.allele=tmp02$mod.locus,
               interact.beta=rowMeans(tmp02[, c('beta', 'beta.replication')], na.rm=T),
               interact.sebeta=rowMeans(tmp02[, c('SEbeta', 'SEbeta.replication')], na.rm=T))
    
  }) %>% do.call(rbind, .)
  out[!is.na(out$mbeta), ]
  
}

haplo.betacomparison <- map(c('DQB1*03:02', 'B*27:05', 'DRB1*04:01'), function(x, y) { alleleVshaplo(x) }) %>% 
  do.call(rbind, .)
haplo.betacomparison$tag[grepl('DM|DIAB|Diab|INSULIN|E4', haplo.betacomparison$pheno)] <- 'Diabetes'
haplo.betacomparison$tag[grepl('DM|DIAB|Diab', haplo.betacomparison$tag)] <- 'Diabetes'
haplo.betacomparison$tag[haplo.betacomparison$pheno=='AUTOIMMUNE'] <- 'Autoimmune'
haplo.betacomparison$tag[haplo.betacomparison$tag=='M13'] <- 'Rheumatoid arthritis'
haplo.betacomparison$tag[haplo.betacomparison$tag=='K11'] <- 'Diseases of the digestive system'
haplo.betacomparison$tag[haplo.betacomparison$tag=='H7'] <- 'Diseases of the eye'
haplo.betacomparison$tag[haplo.betacomparison$tag=='ILD_CM'] <- 'ILD-related Thyroiditis'
haplo.betacomparison$tag[haplo.betacomparison$tag=='RHEUMA'] <- 'Rheumatoid arthritis'
haplo.betacomparison$tag[haplo.betacomparison$tag=='RHEUMA_CM'] <- 'Rheumatoid arthritis'
haplo.betacomparison$tag[haplo.betacomparison$tag=='ASTHMA_CM'] <- 'Asthma'
haplo.betacomparison$tag[haplo.betacomparison$tag=='GASTRO'] <- 'Diseases of the digestive system'
haplo.betacomparison <- filter(haplo.betacomparison, tag!='RX') %>% filter(., tag!='Z21')
haplo.betacomparison$tag2 <- factor(haplo.betacomparison$tag, 
                                   levels=c(geno.betacomparison$tag, haplo.betacomparison$tag) %>% unique)
haplo.betacomparison$diff <- haplo.betacomparison$interact.beta - haplo.betacomparison$comp.beta

haplo.betacomparison.f <- map(haplo.betacomparison$pheno %>% unique, function(x) {
  out <- filter(haplo.betacomparison, pheno==x)
  arrange(out, 1/abs(diff))[1, c(1:2, 4,3, 5:7, 10:11, 15:16, 18)]
}) %>% do.call(rbind, .) %>% arrange(., 1/(abs(diff))) %>% 
  filter(., pheno %in% c('RHEUMA_SERONEG', 'T1D', 'RHEUMA_SEROPOS_STRICT', 'K11_REIMB_202', 'ASTHMA_COMORB' ,'RHEUMA_OTHER_WIDE')) %>% 
  arrange(., diff)
haplo.betacomparison.f$pheno <- paste0('(', haplo.betacomparison.f$pheno ,')')
haplo.betacomparison.f <- unite(haplo.betacomparison.f, fullname, longname,  pheno, sep=' ')
haplo.betacomparison.f <- unite(haplo.betacomparison.f, HLAs, risk.gene, locus, sep=' - ')
haplo.betacomparison.f$fullname <- factor(haplo.betacomparison.f$fullname, levels=haplo.betacomparison.f$fullname %>% unique)

