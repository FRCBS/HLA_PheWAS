

library(tidyverse)
library(data.table)
library(SPAtest)


# ----------------------------------------------
# functions
# ----------------------------------------------

# dosage format
makeDosage <- function(x, gene) {
  # x is HLA imputation output file 
  # gene is HLA gene
  alleles.uniq <- na.omit(unique(unlist(x[, c('allele1', 'allele2')])))
  out <- lapply(alleles.uniq, function(i) {
    apply(x, 1, function(x) as.numeric(x==i)) %>% t %>% rowSums
  }) %>% do.call(cbind, .)
  colnames(out) <- paste0(gene, '*', alleles.uniq)
  data.frame(sample.id=x[, 1], out)

}


# ----------------------------------------------
# read data
# ----------------------------------------------

# hla data
hla.a <- fread('/home/ivm/Documents/R5/results/R5_A_imputed.tsv', data.table=F)
hla.a.dosage <- makeDosage(hla.a, 'A')

hla.b <- fread('/home/ivm/Documents/R5/results/R5_B_imputed.tsv', data.table=F)
hla.b.dosage <- makeDosage(hla.b, 'B')

hla.c <- fread('/home/ivm/Documents/R5/results/R5_C_imputed.tsv', data.table=F)
hla.c.dosage <- makeDosage(hla.c, 'C')

hla.dpb1 <- fread('/home/ivm/Documents/R5/results/R5_DPB1_imputed.tsv', data.table=F)
hla.dpb1.dosage <- makeDosage(hla.dpb1, 'DPB1')

hla.dqa1 <- fread('/home/ivm/Documents/R5/results/R5_DQA1_imputed.tsv', data.table=F)
hla.dqa1.dosage <- makeDosage(hla.dqa1, 'DQA1')

hla.dqb1 <- fread('/home/ivm/Documents/R5/results/R5_DQB1_imputed.tsv', data.table=F)
hla.dqb1.dosage <- makeDosage(hla.dqb1, 'DQB1')

hla.drb1 <- fread('/home/ivm/Documents/R5/results/R5_DRB1_imputed.tsv', data.table=F)
hla.drb1.dosage <- makeDosage(hla.drb1, 'DRB1')

hla.dosage <- Reduce(inner_join, list(hla.a.dosage, hla.b.dosage, hla.c.dosage, hla.dpb1.dosage, hla.dqa1.dosage, hla.dqb1.dosage, hla.drb1.dosage))


# covariates
covs <- fread('./data/R5sansR3_covariates.tsv', data.table=F)
colnames(covs)[1] <- 'sample.id'

# phenotypes
phenos <- fread('./data/R5_pheno_pruned_sansR3.tsv', data.table=F)
phenos <- select(phenos, -c(2:4))
colnames(phenos)[1] <- 'sample.id'



# ----------------------------------------------
# run hla assoc phewas
# ----------------------------------------------

phewas <- lapply(2:ncol(phenos), function(i) {

	print(i)

	# input for spatest
	tmp.cov    <- inner_join(na.omit(phenos[, c(1, i)]), covs, by='sample.id')
	tmp.pheno  <- tmp.cov[, 2]
	tmp.cov    <- tmp.cov[, -2]
	tmp.dosage <- left_join(tmp.cov[, 1:2], hla.dosage)
	tmp.dosage <- tmp.dosage[, -c(1:2)]
	tmp.cov    <- select(tmp.cov, -sample.id) 

	# run test
	tmp <- ScoreTest_SPA(as.matrix(t(tmp.dosage)), tmp.pheno, tmp.cov, beta.out=T, beta.Cutoff=0.05)
	tmp <- do.call(cbind, tmp)
	data.frame(locus=colnames(tmp.dosage), tmp)

})
names(phewas) <- colnames(phenos)[2:ncol(phenos)]

# save result list
save(phewas, file='./results/R5_HLA_PheWAS.RData')



# ----------------------------------------------
# arrange results
# ----------------------------------------------

phewas.all <- map2(phewas, colnames(phenos)[-1], function(x, y) {
	tmp <- data.frame(pheno=y, locus=as.character(x[, 1]), x[, c(2,5,6)], stringsAsFactors=F)
	tmp <- arrange(tmp, p.value)
}) %>% do.call(rbind, .)

phewas.all$locus <- data.frame(str_split_fixed(phewas.all$locus, '\\.', 3)) %>% 
	unite(., TMP, X1,X2, sep='*') %>% unite(., Locus, TMP,X3, sep=':')

phewas.all$p.value <- unlist(phewas.all$p.value)
phewas.all$beta    <- unlist(phewas.all$beta)
phewas.all$SEbeta  <- unlist(phewas.all$SEbeta)
phewas.all$pheno   <- as.character(phewas.all$pheno)
phewas.all$locus   <- phewas.all$locus[, 1]

# write out as a table
write.table(phewas.all, '/home/ivm/Documents/R5/results/R5sansR3_HLA_PheWAS_all_26Mar20.tsv', sep='\t', row.names=F)



