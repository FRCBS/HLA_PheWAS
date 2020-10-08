

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
write.table(haplo.results.z, './results/R3_haplotype_betadiff_replicated_09June20.tsv', sep='\t', row.names=F)

# join with association anno results
# haplo.results.z.anno <- left_join(unite(haplo.results.z[, c(2,4:5,29:31)], phenoallele, pheno, risk.locus, mod.locus, sep='_'),
#                                   unite(haplo.results.anno, phenoallele, pheno, risk.locus, mod.locus, sep='_', remove=F), 
#                                   by='phenoallele') %>% dplyr::select(., -phenoallele) 
# haplo.results.z.anno <- haplo.results.z.anno[, c(4:12, 1:2, 13:23, 3, 24:31)]
# haplo.results.z.anno.replicated <- haplo.results.z.anno[
#   ifelse(haplo.results.z.anno$replication.beta.diff.p.value<0.01 | 
#            is.na(haplo.results.z.anno$replication.beta.diff.p.value), T, F), ]
# haplo.results.z.anno.replicated <- haplo.results.z.anno.replicated %>% add_count(pheno) %>% 
#   filter(n>1) %>% dplyr::select(-n) %>% data.frame
# write 
haplo.results.z.list <- tapply(1:nrow(haplo.results.z), haplo.results.z$pheno, function(x) haplo.results.z[x, ]) 
write.table(lapply(haplo.results.z.list, function(x) 
  rbind(data.frame(x), rep(' ', ncol(haplo.results.z.anno.replicated)))) %>% 
              do.call(rbind, .), './results/tables/PheWAS_haplotype_betadiffs_replication_09June20.tsv', sep='\t', row.names=F)


       
## --------------------------------------------------------
## beta distribution plot
## --------------------------------------------------------

# arrange haplo results by modulation effect magnitude
haplo.results.z.diffs <- (haplo.results.z$pheno %>% table %>% names) %>% map(function(x) {
  tt        <- filter(haplo.results.z, pheno==x, risk.locus.mean.imp.pprob>0.6, risk.locus.mean.imp.pprob.replication>0.6)
  ind.min   <- which.min(tt$beta)
  ind.max   <- which.max(tt$beta)
  beta.min  <- tt[ind.min, c('beta', 'beta.replication')] %>% rowMeans
  beta.max  <- tt[ind.max, c('beta', 'beta.replication')] %>% rowMeans
  p.rep <- tt[ind.min, 'replication.beta.diff.p.value'] 
  data.frame(pheno=x, betadiff=abs(beta.min-beta.max), p.rep=p.rep)
}) %>% do.call(rbind, .)
haplo.results.z.diffs <- arrange(haplo.results.z.diffs, desc(betadiff))
haplo.results.z.diffs <- filter(haplo.results.z.diffs, p.rep<0.01)

# plot all

# make a selection of representative phenos
p.haplo.betas <- ggarrange(
  plotlist=map(c('T1D', 'H7_RETINOPATHYDIAB', 'INSULINS', 
                 'RHEUMA_SEROPOS_STRICT', 'RHEUMA_SERONEG', 'RHEUMA_OTHER_WIDE', 
                 'M13_POLYARTHROPATHIES', 'IBD_DRUGS', 'RX_CROHN_2NDLINE'), function(x) {
    num  <- match(x, c('T1D', 'H7_RETINOPATHYDIAB', 'INSULINS', 
                       'RHEUMA_SEROPOS_STRICT', 'RHEUMA_SERONEG', 'RHEUMA_OTHER_WIDE', 
                       'M13_POLYARTHROPATHIES', 'IBD_DRUGS', 'RX_CROHN_2NDLINE'))
    tmp  <- haplo.results.z[haplo.results.z$pheno==x, ]
    tmp.allele <- tmp[1, 'risk.locus'] # pick risk allele
    tmp.allele.betas <- filter(results.anno.stats.fdr.replicated, 
                               pheno==x, locus==tmp.allele)[1, c('beta', 'beta.replication')]
    tmp$bmean <- rowMeans(dplyr::select(tmp, beta, beta.replication))
    tmp  <- tmp %>% arrange(desc(bmean))
    tmp2 <- data.frame(tmp[, c(2,4:5,19:20)], group='replication')
    tmp  <- data.frame(tmp[, c(2,4:5,8:9)], group='discovery')
    colnames(tmp2) <- colnames(tmp)
    tmp  <- rbind(tmp, tmp2)
    tmp  <- tmp[nrow(tmp):1, ]
    tmp$mod.locus <- factor(tmp$mod.locus, levels=unique(tmp$mod.locus))
    ggplot(tmp, aes(beta, mod.locus, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
      geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.6, fatten=6) +
      ylab(paste0(tmp$risk.locus[1], ' -')) + 
      xlab(ifelse(num==8, 'beta', '')) +
      geom_vline(xintercept=tmp.allele.betas[1, 1], linetype='dashed', alpha=0.7, size=0.3, color="#F8766D") +
      geom_vline(xintercept=tmp.allele.betas[1, 2], linetype='dashed', alpha=0.7, size=0.3, color="#00BFC4") +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.text.y=element_text(size=8.5), 
            axis.line.x=element_line(size=.3), axis.line.y=element_line(size=.3), axis.title.y=element_text(color='black'),
            strip.background=element_rect(fill=NA), legend.position='top',
            legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA),
            legend.text=element_text(size=11))
    }), align='hv', heights=c(1, 0.8, 0.8), common.legend=T, nrow=3, ncol=3
) 
#cairo_pdf('./results/haplo_beta_pointrange_selected_09June20.pdf', height=9, width=9.5)
jpeg('./results/figures/Figure5.jpg', height=9, width=9.5, unit='in', res=600)
p.haplo.betas
dev.off()



# diabetic diseases
cairo_pdf('./results/haplo_beta_pointrange_diab_09June20.pdf', height=34, width=9)
ggarrange(
  plotlist=map(haplo.results.z.diffs$pheno[grepl('DIAB|T1D|T2D|^D1_|^D2_|^DM|E4_D1_|E4_D2_|E4_DM', 
                                                 haplo.results.z.diffs$pheno)], function(x) {
    num  <- which(haplo.results.z.diffs$pheno[grepl('DIAB|T1D|T2D|^D1_|^D2_|^DM|E4_D1_|E4_D2_|E4_DM', 
                                                    haplo.results.z.diffs$pheno)]==x)
    tmp  <- haplo.results.z[haplo.results.z$pheno==x, ]
    tmp2 <- data.frame(tmp[, c(2,4:5,19:20)], group='replication')
    tmp  <- data.frame(tmp[, c(2,4:5,8:9)], group='discovery')
    colnames(tmp2) <- colnames(tmp)
    tmp  <- rbind(tmp, tmp2)
    tmp  <- tmp[nrow(tmp):1, ]
    tmp$mod.locus <- factor(tmp$mod.locus, levels=unique(tmp$mod.locus))
    ggplot(tmp, aes(beta, mod.locus, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
      geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
      ylab(paste0(tmp$risk.locus[1], ' -')) + 
      xlab(ifelse(num==8, 'beta', '')) +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.text.y=element_text(size=8.5), 
            axis.line.x=element_line(size=.2), axis.line.y=element_line(size=.2), axis.title.y=element_text(color='black'),
            strip.background=element_rect(fill=NA), legend.position='top',
            legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA))
  }), align='hv', common.legend=T, nrow=11, ncol=3
) 
dev.off()

# non-diabetics
cairo_pdf('./results/haplo_beta_pointrange_notdiab_09June20.pdf', height=12, width=9.5)
ggarrange(
  plotlist=map(haplo.results.z.diffs$pheno[grepl('DIAB|T1D|T2D|^D1_|^D2_|^DM|E4_D1_|E4_D2_|E4_DM|RX|INSULIN|MED|KELA|DRUG', 
                                                 haplo.results.z.diffs$pheno)==F], function(x) {
    num  <- which(haplo.results.z.diffs$pheno[grepl('DIAB|T1D|T2D|^D1_|^D2_|^DM|E4_D1_|E4_D2_|E4_DM|RX|INSULIN|MED|KELA|DRUG', 
                                                    haplo.results.z.diffs$pheno)==F]==x)
    tmp  <- haplo.results.z[haplo.results.z$pheno==x, ]
    tmp2 <- data.frame(tmp[, c(2,4:5,19:20)], group='replication')
    tmp  <- data.frame(tmp[, c(2,4:5,8:9)], group='discovery')
    colnames(tmp2) <- colnames(tmp)
    tmp  <- rbind(tmp, tmp2)
    tmp  <- tmp[nrow(tmp):1, ]
    tmp$mod.locus <- factor(tmp$mod.locus, levels=unique(tmp$mod.locus))
    ggplot(tmp, aes(beta, mod.locus, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
      geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
      ylab(paste0(tmp$risk.locus[1], ' -')) + 
      xlab(ifelse(num==8, 'beta', '')) +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.text.y=element_text(size=8.5), 
            axis.line.x=element_line(size=.2), axis.line.y=element_line(size=.2), axis.title.y=element_text(color='black'),
            strip.background=element_rect(fill=NA), legend.position='top',
            legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA)) 
  }), align='hv', common.legend=T, nrow=4, ncol=3
) 
dev.off()

# drug purchase
cairo_pdf('./results/haplo_beta_pointrange_medication_09June20.pdf', height=16, width=9.5)
ggarrange(
  plotlist=map(haplo.results.z.diffs$pheno[grepl('RX|INSULIN|MED|KELA|DRUG', haplo.results.z.diffs$pheno)], function(x) {
    num  <- which(haplo.results.z.diffs$pheno[grepl('RX|INSULIN|MED|KELA|DRUG', haplo.results.z.diffs$pheno)]==x)
    tmp  <- haplo.results.z[haplo.results.z$pheno==x, ]
    tmp2 <- data.frame(tmp[, c(2,4:5,19:20)], group='replication')
    tmp  <- data.frame(tmp[, c(2,4:5,8:9)], group='discovery')
    colnames(tmp2) <- colnames(tmp)
    tmp  <- rbind(tmp, tmp2)
    tmp  <- tmp[nrow(tmp):1, ]
    tmp$mod.locus <- factor(tmp$mod.locus, levels=unique(tmp$mod.locus))
    ggplot(tmp, aes(beta, mod.locus, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
      geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
      ylab(paste0(tmp$risk.locus[1], ' -')) + 
      xlab(ifelse(num==8, 'beta', '')) +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.text.y=element_text(size=8.5), 
            axis.line.x=element_line(size=.2), axis.line.y=element_line(size=.2), axis.title.y=element_text(color='black'),
            strip.background=element_rect(fill=NA), legend.position='top',
            legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA))
  }), align='hv', common.legend=T, nrow=5, ncol=3
) 
dev.off()






tmp$bmean <- rowMeans(dplyr::select(tmp, beta, beta.replication))
tmp <- tmp %>% arrange(desc(bmean))

## --------------------------------------------------------
## protein alevel allele divergenses
## --------------------------------------------------------

# HLA protein alignments
hla.aligned <- list(readAAStringSet('~/Projects/gamete2/data/aligned/A_prot_aligned.fasta'),
                    readAAStringSet('~/Projects/gamete2/data/aligned/B_prot_aligned.fasta'),
                    readAAStringSet('~/Projects/gamete2/data/aligned/C_prot_aligned.fasta'),
                    readAAStringSet('~/Projects/gamete2/data/aligned/DRB_prot_aligned.fasta'),
                    readAAStringSet('~/Projects/gamete2/data/aligned/DQA1_prot_aligned.fasta'),
                    readAAStringSet('~/Projects/gamete2/data/aligned/DQB1_prot_aligned.fasta'),
                    readAAStringSet('~/Projects/gamete2/data/aligned/DPB1_prot_aligned.fasta'))
names(hla.aligned) <- c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1')                    

# AA divergence metrics, incl. Grantham
divergences <- fread('~/Projects/gamete2/data/AA_divergence_metrics.tsv', data.table=F)

# divergences of min max alleles

haplo.results.z$risk.locus <- gsub('.', '*', haplo.results.z$risk.locus, fixed=T)
haplo.results.z$mod.locus.max <- str_split_fixed(haplo.results.z$mod.locus.max, '\\.', 3) %>% data.frame %>% 
  unite(., allele, 1, 2, sep='*') %>% unite(., allele, 1, 2, sep=':') %>% .[, 1]
haplo.results.z$mod.locus.min <- str_split_fixed(haplo.results.z$mod.locus.min, '\\.', 3) %>% data.frame %>% 
  unite(., allele, 1, 2, sep='*') %>% unite(., allele, 1, 2, sep=':') %>% .[, 1]

haplo.results.z <- mutate(haplo.results.z, grantham=sapply(1:nrow(haplo.results.z), function(i) {
  HLAdivergence(c(haplo.results.z$mod.locus.max[i], haplo.results.z$mod.locus.min[i]),
                divergences, hla.aligned[[ haplo.results.z$mod.locus.max[i] %>% str_split(., '\\*') %>% .[[1]] %>% .[1] ]], 
                'grantham') })
  )
haplo.results.z <- mutate(haplo.results.z, maxminrel=beta.max/beta.min)

plot(haplo.results.z$grantham, haplo.results.z$maxminrel)
haplo.results.z$grantham %>% hist(xlim=c(0, 15))
haplo.results.z$grantham %>% mean

random.grantham <- sapply(1:2000, function(i) {
  loc <- sample(names(hla.aligned), 1)
  HLAdivergence(sample(names(hla.aligned[[loc]]), 2), divergences, hla.aligned[[loc]], 'grantham') })
random.grantham %>% hist(xlim=c(0, 15))
random.grantham %>% mean


####### tests ############


# correlate haplo freq with beta

fer.in <- fread('../HaploDP/data/new_fers', data.table=F, header=F)
fer.in <- apply(fer.in, 1, gsub, pattern='[ABC]|DQB1|DRB1', replacement='') %>% 
  apply(., 1, gsub, pattern='^\\*|;$|,$', replacement='') %>% data.frame %>% 
  unite(., FER, sep='-', 1:5)
fer.ex <- fread('../HaploDP/data/exclude_fer', data.table=F, header=F, sep=' ')
fer.ex <- apply(fer.ex, 1, gsub, pattern='[ABC]|DQB1|DRB1', replacement='') %>% 
  apply(., 1, gsub, pattern='^\\*|;$|,$', replacement='') %>% data.frame %>% 
  unite(., FER, sep='-', 1:5)
tmp <- fread('../HaploDP/data/FER_haplo_DP', data.table=F, na.strings=c('NA', ''))
tmp <- tmp %>% fill(., 'A-B-C-DR-DQ') # fill missing using previous entry
tmp$`A-B-C-DR-DQ` <- gsub('-$', '', tmp$`A-B-C-DR-DQ`)
tmp$FER <- !(tmp$FER %>% is.na)
tmp$FER <- ifelse(tmp$FER==T, 'Yes', 'No')
tmp$FER[tmp$`A-B-C-DR-DQ` %in% fer.in[, 1]] <- 'Yes'
tmp$FER[tmp$`A-B-C-DR-DQ` %in% fer.ex[, 1]] <- 'No'
tmp <- cbind(tmp, dplyr::select(tmp, "A-B-C-DR-DQ")[,1] %>% str_split_fixed(., '-', 5))
colnames(tmp)[10:14] <- c('A', 'B', 'C', 'DRB1', 'DQB1')

tt <- filter(haplo.results, pheno=='AB1_ARTHROPOD')
tt.1 <- str_split_fixed(tt$risk.locus, '\\*', 2)[1, ]
tt.2 <- str_split_fixed(tt$mod.locus, '\\*', 2)
dplyr::select(tmp, tt.1[1])
filter(tmp, C==tt.1[2], DQB1==tt.2[1,2])

# beta plots by disease 

(haplo.results.z$pheno %>% table %>% names)

intersect((haplo.results.anno$pheno %>% table %>% names), (geno.results.unf$pheno %>% table %>% names))
intersect((haplo.results.z$pheno %>% table %>% names), (geno.results.z$pheno %>% table %>% names))

cairo_pdf('./results/haplo_beta_pointrange_08Apr20.pdf', height=7, width=8)
ggarrange(
  plotlist=map(c('T1D', 'DM_RETINOPATHY', 'DIAB_MED',
                 'RHEUMA_SERONEG', 'RX_CROHN_2NDLINE', 'M13_ARTHROPATHIES', 
                 'IBD_COMORB', 'ASTHMA_COMORB', 'H7_IRIDOCYCLITIS'), function(x) {
    tmp <- haplo.results.z[haplo.results.z$pheno==x, ]
    tmp2 <- data.frame(tmp[, c(2,4:5,15:16)], group='R')
    tmp <- data.frame(tmp[, c(2,4:5,8:9)], group='D')
    colnames(tmp2) <- colnames(tmp)
    tmp <- rbind(tmp, tmp2)
    tmp <- tmp[nrow(tmp):1, ]
    tmp$mod.locus <- factor(tmp$mod.locus, levels=unique(tmp$mod.locus))
    ggplot(tmp, aes(mod.locus, beta, ymin=beta-SEbeta, ymax=beta+SEbeta, color=group)) +
      geom_pointrange(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
      xlab(paste0(tmp$risk.locus[1], ' -')) + 
      ylab(ifelse(x=='ASTHMA_COMORB', 'beta', '')) +
      #ylim(c(-1, 3.7)) +
      #scale_color_manual(values=c("#999999", "blue")) +
      coord_flip() +
      #geom_hline(yintercept=0, linetype='dotted', size=0.3) +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.text.y=element_text(size=8.5), 
            axis.line.x=element_line(size=.2), axis.line.y=element_line(size=.2), axis.title.y=element_text(color='black'),
            strip.background=element_rect(fill=NA), legend.position='none')
  }), align='hv', common.legend=F
) 
dev.off()


cairo_pdf('./results/haplo_beta_dist_flip_02Apr20.pdf', height=8, width=8)
ggarrange(
  plotlist=map(c('T1D', 'RHEUMA_SEROPOS_STRICT', 'M13_ANKYLOSPON', 'L12_PSORIASIS', 'K11_COELIAC', 
                 'STILL_ADULT', 'G6_DEMYEL', 'AB1_SEXUAL_TRANSMISSION', 'H7_VITRHAEMORR'), function(x) {
                   tt <- haplo.results.anno[haplo.results.anno$pheno==x, ] %>% arrange(., beta)
                   #tt$beta[tt$beta<0] <- 0.01
                   tt$mod.locus <- factor(tt$mod.locus, levels=unique(tt$mod.locus))
                   tt$SEbeta[tt$beta<0] <- tt$SEbeta[tt$beta<0]*(-1)
                   ggplot(tt, aes(mod.locus, beta)) +
                     geom_bar(stat='identity', fill='grey70') +
                     geom_errorbar(aes(ymin=beta, ymax=beta+SEbeta), size=.2, width=0.4, position=position_dodge(.9), 
                                   colour='grey70') + 
                     xlab(tt$risk.locus[1]) + ylab(ifelse(x=='M13_ANKYLOSPON', '\nbeta', '')) +
                     #ylim(c(0, 3.7)) +
                     scale_y_continuous(limits=c(ifelse(min(tt$beta+tt$SEbeta)>0, 0, min(tt$beta+tt$SEbeta)), 3.5), expand=c(0, 0)) +
                     coord_flip() +
                     facet_wrap(~pheno, scales='free_x') +
                     theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
                           axis.text.x=element_text(size=9), axis.text.y=element_text(size=8.5), 
                           axis.line.x=element_line(size=.2), axis.title.y=element_text(color='black'),
                           strip.background=element_rect(fill=NA))
                 }), align='hv'
) 
dev.off()




tmp <- haplo.results.z[haplo.results.z$pheno=='M13_RHEUMA', ]
tmp2 <- data.frame(tmp[, c(4:5,15:16)], group='R')
tmp <- data.frame(tmp[, c(4:5,8:9)], group='D')
colnames(tmp2) <- colnames(tmp)
tmp <- rbind(tmp, tmp2); rm(tmp2)
tmp <- tmp[nrow(tmp):1, ]
tmp$mod.locus <- factor(tmp$mod.locus, levels=unique(tmp$mod.locus))

ggplot(tmp, aes(mod.locus, beta, ymin=beta-SEbeta, ymax=beta+SEbeta, col=group)) +
  geom_pointrange(position=position_dodge(width=0.5), shape=18, size=0.6, fatten=6) +
  xlab(tmp$risk.locus[1]) + ylab('beta') +
  ylim(c(-1, 3.7)) +
  coord_flip()



##
haplo.results.z <- tapply(1:nrow(haplo.results.anno), haplo.results.anno$pheno, function(x) {
  tmp <- haplo.results.anno[x, ]
  ind.max <- which.max(tmp$beta)
  ind.min <- which.min(tmp$beta)
  out <- Z.test.2(tmp[ind.max, 'beta'], tmp[ind.min, 'beta'], tmp[ind.max, 'SEbeta'], tmp[ind.min, 'SEbeta'])
  data.frame(pheno=tmp$pheno[1], risk.locus=tmp$risk.locus[1], 
             mod.locus.max=tmp$mod.locus[ind.max], beta.max=tmp$beta[ind.max], beta.max.SE=tmp$SEbeta[ind.max],
             mod.locus.min=tmp$mod.locus[ind.min], beta.min=tmp$beta[ind.min], beta.min.SE=tmp$SEbeta[ind.min],
             p.value=out)
}) %>% do.call(rbind, .)

