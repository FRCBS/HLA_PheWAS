

library(data.table)
library(tidyverse)
library(SPAtest)


## HLA class I - class II haplotype PheWAS


## ----------------------------------------
## Data
## ----------------------------------------

# imputation results
load('./results/R5_haplo.tab.RData')
load('./results/R5_HLA_haplo_imputation.RData')

# allele phewas results
phewas <- fread('./results/PheWAS_fdr001_replication.tsv', data.table=F)
phewas <- phewas[!grepl('\\:$', phewas$locus), ] # remove ambiguous

# covariates
covs <- fread('./data/R5_covariates.tsv', data.table=F)
colnames(covs)[1] <- 'sample.id'

# phenos
phenos <- fread('./data/R5_pheno_pruned.tsv', data.table=F)
phenos <- select(phenos, -c(2:4))
colnames(phenos)[1] <- 'sample.id'

# for each pheno, find the top association and the second best from HLA I or II, if exists
top.assocs <- map(unique(phewas$pheno), function(i) {

  tmp <- filter(phewas, pheno==i) # data subset

  # top risk locus
  loc <- filter(tmp, min(p.value), beta>0) %>% top_n(., -1, p.value) %>% .$locus
  if(length(loc)==0) loc <- NA
  # extract class of risk gene
  tmp.class <- str_split_fixed(loc, '\\*', 2)[1] %>% nchar
  loc2 <- if(tmp.class>1) {
    filter(tmp, grepl('^A|^B|^C', locus))
    } else filter(tmp, !grepl('^A|^B|^C', locus))
  loc2 <- loc2[which.min(tmp$p.val), 'locus']
  
  # output: primary risk allele and secondary allele of a different class
  c(loc, loc2)

})
names(top.assocs) <- unique(phewas$pheno)

# filter out phenos without class I -class II allele pair
# output is a list of allele pairs in haplo loci
top.assocs <- top.assocs[unlist(lapply(top.assocs, function(x) !any(is.na(x))))]


## ----------------------------------------
## run hla haplo assoc phewas
## ----------------------------------------

h.phewas <- lapply(1:length(top.assocs), function(i) {
  
  print(i)
    
  # make an allele table
  all.tab <- str_split_fixed(top.assocs[[i]], '\\*', 2)
  
  # extract haplotypes: risk allele with all available secondary alleles
  tmp.haplo <- haplo.tab[ haplo.tab[, all.tab[1, 1]]==all.tab[1, 2], 
                          c('sample.id', all.tab[1, 1], all.tab[2, 1]) ]
  colnames(tmp.haplo)[2:3] <- c('RISK', 'MOD') 

  # spread into a dosage matrix
  tmp.haplo.mat   <- model.matrix(RISK ~ MOD + 0, data=tmp.haplo)

  # filter out rare allele combinations
  tmp.haplo.mat   <- tmp.haplo.mat[, colSums(tmp.haplo.mat)>19]
  tmp.haplo.mat   <- data.frame(sample.id=tmp.haplo$sample.id, tmp.haplo.mat)
  
  # combine data from both alleles into a single dosage value
  tmp.haplo.mat.1 <- tmp.haplo.mat[grepl('_1', tmp.haplo$sample.id), ]
  tmp.haplo.mat.2 <- tmp.haplo.mat[grepl('_2', tmp.haplo$sample.id), ]
  tmp.haplo.mat.1$sample.id <- gsub('_1|_2', '', tmp.haplo.mat.1$sample.id)
  tmp.haplo.mat.2$sample.id <- gsub('_1|_2', '', tmp.haplo.mat.2$sample.id)
  tmp.haplo.mat <- full_join(tmp.haplo.mat.1, tmp.haplo.mat.2, by='sample.id')
  tmp.haplo.mat[is.na(tmp.haplo.mat)] <- 0
  tmp.haplo.mat <- data.frame(sample.id=tmp.haplo.mat$sample.id, 
                              tmp.haplo.mat[, grepl('x', colnames(tmp.haplo.mat))] + 
                                tmp.haplo.mat[, grepl('y', colnames(tmp.haplo.mat))])
  colnames(tmp.haplo.mat) <- gsub('.x', '', colnames(tmp.haplo.mat), fixed=T)
  colnames(tmp.haplo.mat) <- gsub('MOD', paste0(all.tab[2,1], '.'), 
                                  colnames(tmp.haplo.mat), fixed=T)
  
  # input for spatest
  tmp.pheno  <- select(phenos, c('sample.id', names(top.assocs)[i])) 
  tmp.cov    <- inner_join(na.omit(tmp.pheno), covs, by='sample.id')[, -2]
  tmp.pheno  <- inner_join(tmp.pheno, tmp.cov[, 1:2], by='sample.id')[, -3]
  tmp.dosage <- left_join(tmp.pheno, tmp.haplo.mat, by='sample.id')[, -2]
  tmp.dosage[is.na(tmp.dosage)] <- 0
  
  tmp.pheno  <- select(tmp.pheno, -sample.id)
  tmp.cov    <- select(tmp.cov, -sample.id)
  tmp.dosage <- select(tmp.dosage, -sample.id)
  
  # run the association tests over all allele combinations
  tmp <- sapply(colnames(tmp.dosage), function(gg) {
    tmp.cov2    <- data.frame(tmp.cov, select(tmp.dosage, -gg))
    tmp.dosage2 <- select(tmp.dosage, gg, 1)
    tmp <- ScoreTest_SPA(as.matrix(t(tmp.dosage2)), tmp.pheno, tmp.cov2, 
                         beta.out=T, beta.Cutoff=1)
    tmp <- do.call(cbind, tmp)
    tmp <- data.frame(pheno=names(top.assocs)[i], 
                      risk.locus=paste0(all.tab[1,1], '.', all.tab[1,2]),
                      mod.locus=colnames(tmp.dosage2), 
                      tmp)
    tmp[1, ] 
  })
  t(tmp)
  
})
names(h.phewas) <- names(top.assocs)

# save result list
save(h.phewas, file='./results/R5_HLA_haplo_PheWAS.RData')

# reformat columns
h.phewas <- lapply(1:length(h.phewas), function(i) {
  data.frame(apply(h.phewas[[i]][ ,1:3], 2, function(x) as.character(unlist(x))),
             apply(h.phewas[[i]][ ,4:8], 2, function(x) as.numeric(unlist(x))))
})

# write out as a table
h.phewas.tab <- do.call(rbind, h.phewas) %>% na.omit
write.table(h.phewas.tab, './results/R5_haplo_phewas.tsv', 
            quote=F, row.names=F, sep='\t')


## ----------------------------------------
## result stats
## ----------------------------------------

haplo.n <- lapply(1:length(top.assocs), function(i) {

  print(i)
  
 # make an allele table
  all.tab <- str_split_fixed(top.assocs[[i]], '\\*', 2)
  
  # extract haplotypes: risk allele with all available secondary alleles
  tmp.haplo <- haplo.tab[ haplo.tab[, all.tab[1, 1]]==all.tab[1, 2], 
                          c('sample.id', all.tab[1, 1], all.tab[2, 1]) ]
  colnames(tmp.haplo)[2:3] <- c('RISK', 'MOD') 

  # spread into a dosage matrix
  tmp.haplo.mat   <- model.matrix(RISK ~ MOD + 0, data=tmp.haplo)

  # filter out rare allele combinations
  tmp.haplo.mat   <- tmp.haplo.mat[, colSums(tmp.haplo.mat)>19]
  tmp.haplo.mat   <- data.frame(sample.id=tmp.haplo$sample.id, tmp.haplo.mat)
  
  # combine data from both alleles into a single dosage value
  tmp.haplo.mat.1 <- tmp.haplo.mat[grepl('_1', tmp.haplo$sample.id), ]
  tmp.haplo.mat.2 <- tmp.haplo.mat[grepl('_2', tmp.haplo$sample.id), ]
  tmp.haplo.mat.1$sample.id <- gsub('_1|_2', '', tmp.haplo.mat.1$sample.id)
  tmp.haplo.mat.2$sample.id <- gsub('_1|_2', '', tmp.haplo.mat.2$sample.id)
  tmp.haplo.mat <- full_join(tmp.haplo.mat.1, tmp.haplo.mat.2, by='sample.id')
  tmp.haplo.mat[is.na(tmp.haplo.mat)] <- 0
  tmp.haplo.mat <- data.frame(sample.id=tmp.haplo.mat$sample.id, 
                              tmp.haplo.mat[, grepl('x', colnames(tmp.haplo.mat))] + 
                                tmp.haplo.mat[, grepl('y', colnames(tmp.haplo.mat))])
  colnames(tmp.haplo.mat) <- gsub('.x', '', colnames(tmp.haplo.mat), fixed=T)
  colnames(tmp.haplo.mat) <- gsub('MOD', paste0(all.tab[2,1], '.'), 
                                  colnames(tmp.haplo.mat), fixed=T)
  
  # input
  tmp.pheno  <- select(phenos, c('sample.id', names(top.assocs)[i])) 
  tmp.cov    <- inner_join(na.omit(tmp.pheno), covs, by='sample.id')[, -2]
  tmp.pheno  <- inner_join(tmp.pheno, tmp.cov[, 1:2], by='sample.id')[, -3]
  tmp.dosage <- left_join(tmp.pheno, tmp.haplo.mat, by='sample.id')[, -2]
  tmp.dosage[is.na(tmp.dosage)] <- 0
  
  tmp.pheno  <- select(tmp.pheno, -sample.id)
  tmp.cov    <- select(tmp.cov, -sample.id)
  tmp.dosage <- select(tmp.dosage, -sample.id)
  
  tmp.mod  <- colnames(tmp.dosage)

  # risk and mod genes
  tmp.risk.gene <- str_split_fixed(top.assocs[[i]][1], '\\*', 2)[1] 
  tmp.mod.gene  <- str_split_fixed(tmp.mod[1], '\\.', 2)[1]
       
  tmp.hap.r <- haplodat[[tmp.risk.gene]]
  tmp.hap.m <- haplodat[[tmp.mod.gene]]

  # merge both alleles
  tmp.haplo$sample.id <- gsub('_1|_2', '', tmp.haplo$sample.id)
  tmp.haplo <- left_join(tmp.haplo, tmp.hap.r, by='sample.id')
  tmp.haplo <- left_join(tmp.haplo, tmp.hap.m, by='sample.id')

  # imputation post prob data
  tmp1 <- data.frame(tmp.haplo[tmp.haplo$allele1.x==as.character(tmp.haplo$RISK[1]), c('MOD','prob1.x', 'prob1.y')])
  tmp2 <- data.frame(tmp.haplo[tmp.haplo$allele2.x==as.character(tmp.haplo$RISK[1]), c('MOD','prob2.x', 'prob2.y')])
  colnames(tmp1) <- colnames(tmp2) <- c('mod', 'prob_risk', 'prob_mod')
  tmp.probs      <- rbind(tmp1, tmp2)
  tmp.mod        <- unite(str_split_fixed(tmp.mod, '\\.', 3)[, 2:3] %>% data.frame, hh, sep=':')[, 1]
  # mean and sd of PPs of risk and mod alleles
  tmp.probs      <- sapply(tmp.mod, function(x) {
    out <- filter(tmp.probs, mod==x)[, 2:3]
    c(mean(out[, 1]), sd(out[, 1]), mean(out[, 2]), sd(out[, 2]))
  }) %>% t 

  # output table
  data.frame(
    pheno                     = names(top.assocs[i]),
    total.n                   = nrow(na.omit(tmp.pheno)),
    cases.n                   = sum(tmp.pheno),
    haplo.risk.allele         = top.assocs[[i]][1],
    haplo.mod.allele          = colnames(tmp.dosage),
    haplo.risk.allele.mean.pp = tmp.probs[, 1],
    haplo.risk.allele.sd.pp   = tmp.probs[, 2],
    haplo.mod.allele.mean.pp  = tmp.probs[, 3],
    haplo.mod.allele.sd.pp    = tmp.probs[, 4],
    haplo.n                   = colSums(tmp.dosage),
    haplo.n.cases             = colSums(tmp.dosage*tmp.pheno[, 1])
  )
  
}) %>% do.call(rbind, .)

write.table(haplo.n, './results/R5_haplo_stats.tsv', sep='\t', row.names=F)
write.table(filter(haplo.n, haplo.n.cases>4), 
            './results/R5_haplo_stats_n5.tsv', sep='\t', row.names=F)

