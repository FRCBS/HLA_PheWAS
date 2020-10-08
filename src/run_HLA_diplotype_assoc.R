
library(data.table)
library(tidyverse)
library(SPAtest)


## HLA allele diplotype PheWAS


## ----------------------------------------
## read data
## ----------------------------------------

# allele phewas results, pick top allele from each phenotype
phewas <- fread('./results/PheWAS_fdr001_replication.tsv', data.table=F)
# omit small betas
phewas <- group_by(phewas, pheno) %>% top_n(., -1, p.value) %>%
  filter(., abs(beta)>0.3)
phewas <- phewas[!grepl('\\:$', phewas$locus), ] # remove ambiguous

# imputed allele lists
ls.imputed <- list.files('./results', '_imputed.tsv', full.names=T)
imp        <- map(ls.imputed, function(x) fread(x, data.table=F))
names(imp) <- gsub('_imputed.tsv|./results/R5_', '', ls.imputed)

# covariates
covs <- fread('./data/R5_covariates.tsv', data.table=F)
colnames(covs)[1] <- 'sample.id'

# phenos
phenos <- fread('./data/R5_pheno_pruned.tsv', data.table=F)
phenos <- select(phenos, -c(2:4))
colnames(phenos)[1] <- 'sample.id'



## ----------------------------------------
## run hla allele diplotype assoc phewas
## ----------------------------------------

g.phewas <- lapply(phewas$pheno, function(i) {
  
  print(i)
  
  # select a disease and its risk gene and allele
  tmp.locus <- filter(phewas, pheno==i)$locus %>% str_split(., '\\*') %>% .[[1]]
  imp.loci  <- filter(imp[[tmp.locus[1]]], allele1==tmp.locus[2] | allele2==tmp.locus[2])

  # spread allele combinations into a dosage matrix
  imp.mm <- mutate(imp.loci, diplo=paste0(imp.loci$allele1, '-', imp.loci$allele2))
  if(length(unique(imp.mm$diplo))>1) {
    imp.mm <- model.matrix(~ 0 + diplo, imp.mm)
  } else imp.mm <- rep(1, nrow(imp.mm))         
  imp.mm <- data.frame(sample.id=imp.loci$sample.id, imp.mm) 
    
  # input for spatest
  tmp.phe    <- select(phenos, c('sample.id', i)) 
  tmp.cov    <- inner_join(na.omit(tmp.phe), covs, by='sample.id')
  tmp.pheno  <- tmp.cov[, 2]
  tmp.cov    <- tmp.cov[, -2]
  tmp.dosage <- left_join(tmp.cov[, 1:2], imp.mm)
  tmp.dosage <- tmp.dosage[, -c(1:2)]
  tmp.dosage[is.na(tmp.dosage)] <- 0 
  tmp.cov    <- select(tmp.cov, -sample.id)
  
  # run tests over all diplotypes
  tmp <- sapply(colnames(tmp.dosage), function(gg) {
    tmp.cov    <- data.frame(tmp.cov, select(tmp.dosage, -gg))
    tmp.dosage <- select(tmp.dosage, gg, 1)
    tmp <- ScoreTest_SPA(as.matrix(t(tmp.dosage)), tmp.pheno, tmp.cov, beta.out=T, beta.Cutoff=1)
    tmp <- do.call(cbind, tmp)
    tmp <- data.frame(pheno=i, locus=colnames(tmp.dosage), tmp)
    tmp$locus <- gsub('diplo', paste0(tmp.locus[1], '*'), tmp$locus)
    tmp[1, ] 
  })
  t(tmp)
  
})
names(g.phewas) <- phewas$pheno

# save result
save(g.phewas, file='./results/R5_HLA_diplo_PheWAS.RData')




## ----------------------------------------
## results stats
## ----------------------------------------

geno.n <- lapply(phewas$pheno, function(i) {
  
  print(i)
  
  # select a disease and its risk gene and allele
  tmp.locus <- filter(phewas, pheno==i)$locus %>% str_split(., '\\*') %>% .[[1]]
  imp.loci  <- filter(imp[[tmp.locus[1]]], allele1==tmp.locus[2] | allele2==tmp.locus[2])

  # spread allele combinations into a dosage matrix
  imp.mm <- mutate(imp.loci, geno=paste0(imp.loci$allele1, '-', imp.loci$allele2))
  if(length(unique(imp.mm$geno))>1) {
    imp.mm <- model.matrix(~ 0+geno, imp.mm)
  } else imp.mm <- rep(1, nrow(imp.mm))         
  imp.mm <- data.frame(sample.id=imp.loci$sample.id, imp.mm) 

  # extract secondary alleles
  tmp.mod.alleles <- gsub('geno', '', colnames(imp.mm)[-1])
  tmp.mod.alleles <- gsub(gsub(':', '.', tmp.locus[2], fixed=T), 
      '', tmp.mod.alleles, fixed=T)
  tmp.mod.alleles[tmp.mod.alleles=='.'] <- gsub(':', '.', tmp.locus[2], fixed=T)
  tmp.mod.alleles <- gsub('^\\.|\\.?$', '', tmp.mod.alleles)

  # PP probs
  tmp.imp.probs <- imp.loci$prob*imp.mm[, -1]
  tmp.imp.probs <- apply(tmp.imp.probs, 2, function(x) {
    c(mean(x[x!=0], na.rm=T), sd(x[x!=0], na.rm=T))
  }) %>% t

  # arrange output
  tmp.phe    <- select(phenos, c('sample.id', i)) 
  tmp.cov    <- inner_join(na.omit(tmp.phe), covs, by='sample.id')
  tmp.pheno  <- tmp.cov[, 2]
  tmp.cov    <- tmp.cov[, -2]
  tmp.dosage <- left_join(tmp.cov[, 1:2], imp.mm)
  tmp.dosage <- tmp.dosage[, -c(1:2)]
  tmp.dosage[is.na(tmp.dosage)] <- 0 
  tmp.cov    <- select(tmp.cov, -sample.id)

  # output table
  out <- data.frame(pheno          = i,
                    locus          = paste0(tmp.locus, collapse='.'), 
                    mod.allele     = tmp.mod.alleles,
                    imp.pprob.mean = tmp.imp.probs[, 1],
                    imp.pprob.sd   = tmp.imp.probs[, 2],
		    total.n        = nrow(tmp.dosage), 
                    cases.n        = sum(tmp.pheno),
                    geno.alleles   = names(colSums(tmp.dosage)),
                    geno.n         = colSums(tmp.dosage),
                    geno.cases     = colSums(tmp.dosage*tmp.pheno))
  out$geno.alleles <- gsub('geno', '', out$geno.alleles)
  out
  
}) %>% do.call(rbind, .)

write.table(geno.n, './results/R5_geno_stats.tsv', sep='\t', row.names=F)
write.table(filter(geno.n, geno.cases>4), './results/R5_geno_stats_n5.tsv', sep='\t', row.names=F)

