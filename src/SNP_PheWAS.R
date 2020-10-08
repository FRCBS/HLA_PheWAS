
source('./src/functions.R')

theme_set(theme_grey())

## -----------------------------------------------------------------
## Process FinnGen summary stat files
## -----------------------------------------------------------------

# list files
ls.dat <- list.files('./data/file_download_21082020_R5_mhc_sumstats', 'MHC', full.names=T)
# extract chr 6
lapply(ls.dat, function(x) {
  tmp <- fread(x, data.table=F, check.names=T) %>% filter(., X.chrom==6)
  write.table(tmp, x, sep='\t', quote=F, row.names=F)
})


## -----------------------------------------------------------------
## MHC genes, hg38 Manhattan plot
## -----------------------------------------------------------------

# MHC genes
ensembl   <- useDataset("hsapiens_gene_ensembl", mart=useMart(biomart="ENSEMBL_MART_ENSEMBL"))
mhc.genes <- getBM(attributes=c('wikigene_name', 'wikigene_description', 'ensembl_gene_id', 
                                'chromosome_name','start_position', 'end_position') , 
                   filters=c('chromosome_name', 'start', 'end'), 
                   values=list('6', 28510120, 33480577), mart=ensembl) %>% filter(., wikigene_name!='')
mhc.genes <- mhc.genes[!duplicated(mhc.genes$wikigene_name), ]
mhc.genes[, c('start_position', 'end_position')] <- mhc.genes[, c('start_position', 'end_position')]/1e6 
mhc.genes$gene_mean <- apply(mhc.genes[, c('start_position', 'end_position')], 1, mean)
mhc.genes.selected  <- filter(mhc.genes, mhc.genes$wikigene_name %in% 
                                c('HLA-F', 'HLA-G', 'HLA-B', 'HLA-A', 'HLA-E', 'IER3', 'DDR1', 'MUC22', 'HCG27', 'HLA-C', 
                                  'MICA', 'MICB', 'TNF', 'SLC44A4', 'CYP21A2', 'HFE',
                                  'C4A', 'NOTCH4', 'HLA-DRA', 'HLA-DRB1', 'HLA-DQB1', 'TAP2', 'HLA-DPB1') )

# list snp summary stat files
ls.snps <- list.files('./data/sumstats_MHC_finngen_R5_all', 'MHC', full.names=T)
ls.snps.names <- gsub('./data/sumstats_MHC_finngen_R5_all/MHC_finngen_R5_', '', ls.snps, fixed=T)


# selected phenos
target.phenos.ind1 <- ls.snps.names %in% c('AB1_BACTINF_NOS', 'AB1_PARASITIC_NOS', 'AB1_VIRAL_OTHER_INTEST_INFECTIONS')
target.phenos.anno1 <- data.frame(label=c('sole association with DQA1*03:01', 
                                          'sole association with B*27:05', 
                                          'sole association with DQA1*03:01'),
                                  pheno=c('AB1_BACTINF_NOS', 'AB1_PARASITIC_NOS', 'AB1_VIRAL_OTHER_INTEST_INFECTIONS'))
p.snps01 <- multiManhattan(target.phenos.ind1, target.phenos.anno1)

target.phenos.ind2 <- ls.snps.names %in% c('RHEUMA_SEROPOS_STRICT' ,'RHEUMA_SERONEG')
target.phenos.anno2 <- data.frame(label=c('primary association with DRB1*04:01, modified by B*44:02', 
                                          'primary association with B*27:05, modified by DRB1*08:01'),
                                  pheno=c('RHEUMA_SEROPOS_STRICT' ,'RHEUMA_SERONEG'))
p.snps02 <- multiManhattan(target.phenos.ind2, target.phenos.anno2)



