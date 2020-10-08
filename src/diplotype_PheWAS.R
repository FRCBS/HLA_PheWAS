

## Analysis of HLA genotype PheWAS results

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
# geno.results.anno.list <- tapply(1:nrow(geno.results.anno), geno.results.anno$pheno, function(x) geno.results.anno[x, ]) 
# write.table(lapply(geno.results.anno.list, function(x) rbind(data.frame(x), rep(' ', ncol(geno.results.anno)))) %>% 
#               do.call(rbind, .), './results/R3_genotype_PheWAS_full_fdr001_08May20.tsv', sep='\t', row.names=F)
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
## beta distribution plot
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


# geno.results.z.diffs <- (geno.results.z$pheno %>% table %>% names) %>% map(function(x) {
#   tt       <- filter(geno.results.z, pheno==x, imp.pprob.mean>0.6, imp.pprob.mean.replication>0.6)
#   ind.min  <- which.min(tt$adjusted.beta.diff.p.value)
#   ind.max  <- which.max(tt$beta)
#   beta.min <- tt[ind.min, c('beta', 'beta.replication')] %>% rowMeans
#   beta.max <- tt[ind.max, c('beta', 'beta.replication')] %>% rowMeans
#   p.rep    <- tt[ind.min, 'replication.beta.diff.p.value'] 
#   data.frame(pheno=x, betadiff=abs(beta.min-beta.max), p.rep=p.rep)
# }) %>% do.call(rbind, .)
# geno.results.z.diffs <- arrange(geno.results.z.diffs, desc(betadiff))
# geno.results.z.diffs <- filter(geno.results.z.diffs, p.rep<0.01)

## literature search for protective alleles

# mechanisms of protection
# mangalam:HLA class II molecules influence susceptibility versus protection in inflammatory diseases by determining the cytokine profile 
# furukawa: shared motif in peptide and HLA molecule
# Ooi: Treg

# t1d https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3233362/
# dqb1 0603? ei uusi
# https://care.diabetesjournals.org/content/27/3/676 
# https://pubmed.ncbi.nlm.nih.gov/11061540/
# our analysis captured the previously described allele effects

# celiac
# dqb1:03:01? ei löydy selkeää infoa protectiivisuudesta

# seropositive rheuma
# drb1 0801= ei löydy selkeää infoa

# MS
# https://pubmed.ncbi.nlm.nih.gov/20207784/ 0701 is protective
# 0101 protective

# lichen planus
# dqb1 0301 ?

# psoriatic arthropaties
# https://pubmed.ncbi.nlm.nih.gov/17166285/ C 0701 protective

# RX_A07EC01 Sulfasalazine
# hla-b51:01

# make a non-redundant selection of representative phenos
# selected endpoints
# geno.results.z.diffs$pheno[
#   match(c('T1D', 'DM_KIDNEYFAIL', 'INSULINS', 
#           'K11_COELIAC','RHEUMA_SEROPOS_STRICT','G6_DEMYEL', 
#           'L12_PSORI_ARTHRO', 'L12_LICHENPLANUS', 'RX_A07EC01'), 
#         geno.results.z.diffs$pheno)]
p.geno.betas <- ggarrange(
  plotlist=map(c('T1D', 'DM_KIDNEYFAIL', 'INSULINS', 
                 'K11_COELIAC','RHEUMA_SEROPOS_STRICT','G6_DEMYEL', 
                 'L12_PSORI_ARTHRO', 'L12_LICHENPLANUS', 'RX_A07EC01'), function(x) {
            tmp <- geno.results.z[geno.results.z$pheno==x, ] 
            tmp.allele <- tmp[1, c('locus', 'risk.allele')] %>% paste(., collapse='*') # pick risk allele
            tmp.allele.betas <- filter(results.anno.stats.fdr.replicated, 
                                       pheno==x, locus==tmp.allele)[1, c('beta', 'beta.replication')]
            tmp$bmean <- rowMeans(dplyr::select(tmp, beta, beta.replication))
            tmp <- tmp %>% arrange(desc(bmean))
            tmp2 <- data.frame(tmp[, c(2,4:6,18:19)], group='replication')
            tmp  <- data.frame(tmp[, c(2,4:6,9:10)], group='discovery')
            colnames(tmp2) <- colnames(tmp)
            tmp <- rbind(tmp, tmp2)
            tmp <- tmp[nrow(tmp):1, ]
            tmp$mod.allele <- factor(tmp$mod.allele, levels=unique(tmp$mod.allele))
            risk.col <- rep('black', nrow(tmp))
            risk.col[tmp$mod.allele==tmp$risk.allele] <- 'red'
            ggplot(tmp, aes(beta, mod.allele, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
              geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.6, fatten=6) +
              ylab(paste0(tmp$locus[1], '*', tmp$risk.allele[1])) + 
              xlab(ifelse(x=='L12_LICHENPLANUS', '\nbeta', '')) +
              geom_vline(xintercept=tmp.allele.betas[1, 1], linetype='dashed', alpha=0.7, size=0.3, color="#F8766D") +
              geom_vline(xintercept=tmp.allele.betas[1, 2], linetype='dashed', alpha=0.7, size=0.3, color="#00BFC4") +
              facet_wrap(~pheno, scales='free') +
              theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
                    axis.text.x=element_text(size=8.5), axis.line.y=element_line(size=.3),
                    axis.text.y=element_text(size=8.5, color=risk.col), axis.line.x=element_line(size=.3),
                    strip.background=element_rect(fill=NA), legend.position='top',
                    legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA),
                    legend.text=element_text(size=11))
          }), align='hv', ncol=3, nrow=3, heights=c(1, 0.9, 0.6), common.legend=T
) 
#cairo_pdf('./results/geno_beta_pointrange_selected_12Jun20.pdf', height=9, width=9)
jpeg('./results/figures/Figure4.jpg', height=10, width=9, unit='in', res=600)
p.geno.betas
dev.off()




# diabetic
cairo_pdf('./results/geno_beta_pointrange_diab_12Jun20.pdf', height=55, width=9)
ggarrange(
  plotlist=map(geno.results.z.diffs$pheno[grepl('DIAB|T1D|T2D|^D1_|^D2_|^DM|E4_D1_|E4_D2_|E4_DM', 
                                                geno.results.z.diffs$pheno)], function(x) {
    tmp <- geno.results.z[geno.results.z$pheno==x, ]
    tmp$bmean <- rowMeans(dplyr::select(tmp, beta, beta.replication))
    tmp <- tmp %>% arrange(desc(bmean))
    tmp2 <- data.frame(tmp[, c(2,4:6,18:19)], group='replication')
    tmp  <- data.frame(tmp[, c(2,4:6,9:10)], group='discovery')
    colnames(tmp2) <- colnames(tmp)
    tmp <- rbind(tmp, tmp2)
    tmp <- tmp[nrow(tmp):1, ]
    tmp$mod.allele <- factor(tmp$mod.allele, levels=unique(tmp$mod.allele))
    risk.col <- rep('black', nrow(tmp))
    risk.col[tmp$mod.allele==tmp$risk.allele] <- 'red'
    ggplot(tmp, aes(beta, mod.allele, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
      geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
      ylab(paste0(tmp$locus[1], '*', tmp$risk.allele[1])) + xlab('') +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.line.y=element_line(size=.2),
            axis.text.y=element_text(size=8.5, color=risk.col), axis.line.x=element_line(size=.2),
            strip.background=element_rect(fill=NA), legend.position='top',
            legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA))
  }), align='hv', ncol=3, nrow=18, common.legend=T
) 
dev.off()

# non-diabetic
cairo_pdf('./results/geno_beta_pointrange_notdiab_12Jun20.pdf', height=31, width=9)
ggarrange(
  plotlist=map(geno.results.z.diffs$pheno[grepl('DIAB|T1D|T2D|^D1_|^D2_|^DM|E4_D1_|E4_D2_|E4_DM|RX|INSULIN|MED|KELA|DRUG', 
                                                geno.results.z.diffs$pheno)==F], function(x) {
    tmp <- geno.results.z[geno.results.z$pheno==x, ] 
    tmp$bmean <- rowMeans(dplyr::select(tmp, beta, beta.replication))
    tmp <- tmp %>% arrange(desc(bmean))
    tmp2 <- data.frame(tmp[, c(2,4:6,18:19)], group='replication')
    tmp  <- data.frame(tmp[, c(2,4:6,9:10)], group='discovery')
    colnames(tmp2) <- colnames(tmp)
    tmp <- rbind(tmp, tmp2)
    tmp <- tmp[nrow(tmp):1, ]
    tmp$mod.allele <- factor(tmp$mod.allele, levels=unique(tmp$mod.allele))
    risk.col <- rep('black', nrow(tmp))
    risk.col[tmp$mod.allele==tmp$risk.allele[1]] <- 'red'
    ggplot(tmp, aes(beta, mod.allele, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
      geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
      ylab(paste0(tmp$locus[1], '*', tmp$risk.allele[1])) + xlab('') +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.line.y=element_line(size=.2),
            axis.text.y=element_text(size=8.5, color=risk.col), axis.line.x=element_line(size=.2),
            strip.background=element_rect(fill=NA), legend.position='top',
            legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA))
  }), align='hv', ncol=3, nrow=10, common.legend=T
) 
dev.off()

# drug purchase
cairo_pdf('./results/geno_beta_pointrange_medication_12Jun20.pdf', height=10, width=9)
ggarrange(
  plotlist=map(geno.results.z.diffs$pheno[grepl('RX|INSULIN|MED|KELA|DRUG', geno.results.z.diffs$pheno)], function(x) {
    tmp <- geno.results.z[geno.results.z$pheno==x, ] 
    tmp$bmean <- rowMeans(dplyr::select(tmp, beta, beta.replication))
    tmp <- tmp %>% arrange(desc(bmean))
    tmp2 <- data.frame(tmp[, c(2,4:6,18:19)], group='replication')
    tmp  <- data.frame(tmp[, c(2,4:6,9:10)], group='discovery')
    colnames(tmp2) <- colnames(tmp)
    tmp <- rbind(tmp, tmp2)
    tmp <- tmp[nrow(tmp):1, ]
    tmp$mod.allele <- factor(tmp$mod.allele, levels=unique(tmp$mod.allele))
    risk.col <- rep('black', nrow(tmp))
    risk.col[tmp$mod.allele==tmp$risk.allele[1]] <- 'red'
    ggplot(tmp, aes(beta, mod.allele, xmin=beta-SEbeta, xmax=beta+SEbeta, color=group)) +
      geom_pointrangeh(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
      ylab(paste0(tmp$locus[1], '*', tmp$risk.allele[1])) + xlab('') +
      facet_wrap(~pheno, scales='free') +
      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
            axis.text.x=element_text(size=8.5), axis.line.y=element_line(size=.2),
            axis.text.y=element_text(size=8.5, color=risk.col), axis.line.x=element_line(size=.2),
            strip.background=element_rect(fill=NA), legend.position='top',
            legend.title=element_blank(), legend.key=element_rect(fill=NA, color=NA))
  }), align='hv', ncol=3, nrow=3, common.legend=T
) 
dev.off()




# cairo_pdf('./results/geno_beta_pointrange_08Apr29.pdf', height=7, width=8)
# ggarrange(
#   plotlist=map(c('T1D', 'K11_COELIAC', 
#                  'RHEUMA_SEROPOS_WIDE', 'L12_BULLOUS', 'GIANT_CELL_TEMP_ARTERITIS', 
#                  'M13_WEGENER', 'L12_PSORIASIS', 'G6_MS', 'CD2_FOLLICULAR_LYMPHOMA'), function(x) {
#                    tmp <- geno.results.z[geno.results.z$pheno==x, ] %>% arrange(desc(beta))
#                    tmp2 <- data.frame(tmp[, c(2,4:6,16:17)], group='R')
#                    tmp <- data.frame(tmp[, c(2,4:6,9:10)], group='D')
#                    colnames(tmp2) <- colnames(tmp)
#                    tmp <- rbind(tmp, tmp2)
#                    tmp <- tmp[nrow(tmp):1, ]
#                    tmp$mod.allele <- factor(tmp$mod.allele, levels=unique(tmp$mod.allele))
#                    risk.col <- rep('black', nrow(tmp))
#                    risk.col[tmp$mod.allele==tmp$risk.allele] <- 'red'
#                    ggplot(tmp, aes(mod.allele, beta, ymin=beta-SEbeta, ymax=beta+SEbeta, color=group)) +
#                      geom_pointrange(position=position_dodge(width=0.55), shape=18, size=0.5, fatten=5) +
#                      xlab(paste0(tmp$locus[1], '*', tmp$risk.allele[1])) + ylab(ifelse(x=='G6_MS', 'beta', '')) +
#                      #ylim(c(-0.5, 3.4)) +
#                      coord_flip() +
#                      #geom_hline(yintercept=0, linetype='dotted', size=0.3) +
#                      facet_wrap(~pheno, scales='free_x') +
#                      theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
#                            axis.text.x=element_text(size=8.5), axis.line.y=element_line(size=.2),
#                            axis.text.y=element_text(size=8.5, color=risk.col), axis.line.x=element_line(size=.2),
#                            strip.background=element_rect(fill=NA), legend.position='none')
# 
#                  }), align='hv', heights=c(1, 0.5, 0.5)
# ) 
# dev.off()



## --------------------------------------------------------
## single alleles vs. genotype
## --------------------------------------------------------

geno.allele.betas <- map(unique(geno.results.z$pheno), function(x) {
  #print(x)
  tmp <- filter(geno.results.z, pheno==x) %>% mutate(., new.locus=paste(locus, mod.allele, sep='*'))
  tmp <- inner_join(filter(results.anno, pheno==x), tmp, by=c('locus'='new.locus'))
  if(nrow(tmp)>0 & !all(is.na(tmp$beta.x))) {
    tmp1 <- tmp[!is.na(tmp$beta.x), c('pheno.x', 'locus', 'beta.x', 'beta.y')]
    tmp1$dr <- 'D'
    tmp2 <- tmp[!is.na(tmp$beta.x), c('pheno.x', 'locus', 'beta.replication.x', 'beta.replication.y')]
    tmp2$dr <- 'R'
    colnames(tmp2) <- colnames(tmp1)
    rbind(tmp1, tmp2)
  } else NA
})
geno.allele.betas <- geno.allele.betas[!is.na(geno.allele.betas)]
geno.allele.betas <- geno.allele.betas[sapply(geno.allele.betas, nrow)>2]
geno.allele.betas <- do.call(rbind, geno.allele.betas)
geno.allele.betas$group <- rep(NA, nrow(geno.allele.betas))
geno.allele.betas$group[grepl('T1D$', geno.allele.betas$pheno.x)] <- 'T1D'
geno.allele.betas$group[grepl('DM|DIAB', geno.allele.betas$pheno.x)] <- 'Diabetes'
geno.allele.betas$group[grepl('COELIAC', geno.allele.betas$pheno.x)] <- 'CD'
geno.allele.betas$group[grepl('RETINA|IRIDO|VITR|MACULO', geno.allele.betas$pheno.x)] <- 'Eye'
geno.allele.betas$group[grepl('RHEUMA', geno.allele.betas$pheno.x)] <- 'Rheuma'
geno.allele.betas$group[grepl('ARTH', geno.allele.betas$pheno.x)] <- 'Arthritis'
geno.allele.betas$group[is.na(geno.allele.betas$group)] <- 'Other'

cairo_pdf('./results/allele_vs_geno_betas_12May20.pdf', height=16, width=16)
ggplot(geno.allele.betas, aes(beta.x, beta.y)) +
  geom_point(aes(color=dr, shape=dr), alpha=0.9, size=2) +
  geom_smooth(method='lm', formula=y~x, se=F, color='black', size=0.2, linetype='dashed') + 
  xlab('allele beta') + ylab('genotype beta') +
  facet_wrap(~pheno.x, scales='free') +
  theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
        axis.text.x=element_text(size=8.5), axis.line.y=element_line(size=.2),
        axis.text.y=element_text(size=8.5), axis.line.x=element_line(size=.2),
        strip.background=element_rect(fill=NA), legend.position='top')
dev.off()


cairo_pdf('./results/allele_vs_geno_betas_15Apr20.pdf', height=5, width=7)
ggplot(filter(geno.allele.betas, group!='Other'), aes(beta.x, beta.y, color=group)) +
  geom_point(alpha=0.8) +
  geom_smooth(method='lm', se=F, color='black', size=0.2, linetype='dashed') + 
  xlab('allele beta') + ylab('genotype beta') +
  facet_wrap(~group, scales='free') +
  theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
        axis.text.x=element_text(size=8.5), axis.line.y=element_line(size=.2),
        axis.text.y=element_text(size=8.5), axis.line.x=element_line(size=.2),
        strip.background=element_rect(fill=NA), legend.position='none')
dev.off()



tmp <- filter(geno.results.z, pheno=='T1D') %>% mutate(., new.locus=paste(locus, mod.allele, sep='*'))
tmp <- inner_join(filter(results, pheno=='T1D'), tmp, by=c('locus'='new.locus'))
plot(tmp$beta.x, tmp$beta.y)

tmp <- filter(geno.results.z, pheno=='K11_COELIAC') %>% mutate(., new.locus=paste(locus, mod.allele, sep='*'))
tmp <- inner_join(filter(results, pheno=='K11_COELIAC'), tmp, by=c('locus'='new.locus'))
plot(tmp$beta.x, tmp$beta.y)
# DQB1*06:03

tmp <- filter(geno.results.z, pheno=='RHEUMA_SEROPOS_WIDE') %>% mutate(., new.locus=paste(locus, mod.allele, sep='*'))
tmp <- inner_join(filter(results, pheno=='RHEUMA_SEROPOS_WIDE'), tmp, by=c('locus'='new.locus'))
plot(tmp$beta.x, tmp$beta.y)

tmp <- filter(geno.results.z, pheno=='GIANT_CELL_TEMP_ARTERITIS') %>% mutate(., new.locus=paste(locus, mod.allele, sep='*'))
tmp <- inner_join(filter(results, pheno=='GIANT_CELL_TEMP_ARTERITIS'), tmp, by=c('locus'='new.locus'))
plot(tmp$beta.x, tmp$beta.y)



## --------------------------------------------------------
## allele interaction table
## --------------------------------------------------------

geno.results.z$mod.allele.max %>% table



## --------------------------------------------------------
## protein level allele divergenes
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

geno.results.z$mod.allele.max <- unite(geno.results.z, allele, locus, mod.allele.max, sep='*') %>% .$allele
geno.results.z$mod.allele.min <- unite(geno.results.z, allele, locus, mod.allele.min, sep='*') %>% .$allele
geno.results.z <- mutate(geno.results.z, grantham=sapply(1:nrow(geno.results.z), function(i) {
  HLAdivergence(c(geno.results.z$mod.allele.max[i], geno.results.z$mod.allele.min[i]),
                divergences, hla.aligned[[ geno.results.z$locus[i] ]], 'grantham') })
)
geno.results.z <- mutate(geno.results.z, maxminrel=beta.max/beta.min)

plot(geno.results.z$grantham, geno.results.z$maxminrel)
geno.results.z$grantham %>% hist(xlim=c(0, 15))
geno.results.z$grantham %>% mean

# generate random divergences
random.grantham <- sapply(1:10000, function(i) {
  loc <- sample(names(hla.aligned), 1)
  HLAdivergence(sample(names(hla.aligned[[loc]]), 2), divergences, hla.aligned[[loc]], 'grantham') })
random.grantham %>% hist(xlim=c(0, 15))
random.grantham %>% mean

boxplot(list(random.grantham, geno.results.z$grantham))
t.test(random.grantham, geno.results.z$grantham)
plot(geno.results.z$beta.max, geno.results.z$beta.min)

# plot densities of grantham between min and max alleles
allele.comp <- rbind(data.frame(grantham=geno.results.z$grantham, group='genotype'),
                     data.frame(grantham=hap.results.z$grantham, group='haplotype'),
                     data.frame(grantham=random.grantham, group='random'))

p01 <- ggplot(allele.comp, aes(x=grantham, group=group)) +
  geom_density(aes(fill=group, color=group), alpha=0.5) +
  xlab("Grantham's distance between modulatory alleles") +
  theme(panel.background=element_blank(), axis.line=element_blank(),#element_line(colour="white", size=.2),
        panel.border=element_rect(fill=NA, size=.2), legend.position=c(.85, .85), legend.title=element_blank(),
        axis.text.x=element_text(angle=0, hjust=1))

cairo_pdf('./results/grantham_dists_16Mar20.pdf', height=4, width=5)
p01
dev.off()

# grantham between mod. allele and risk allele  
allele.div <- filter(geno.results.anno, pheno %in% geno.results.z$pheno)
allele.div$risk.allele <- unite(allele.div, tt, locus, risk.allele, sep='*')$tt
allele.div$mod.allele  <- unite(allele.div, tt, locus, mod.allele, sep='*')$tt
allele.div <- mutate(allele.div, grantham=sapply(1:nrow(allele.div), function(i) {
  HLAdivergence(c(allele.div$risk.allele[i], allele.div$mod.allele[i]),
                divergences, hla.aligned[[ allele.div$locus[i] ]], 'grantham') })
) 
# select min and max mod alleles and their grantham distances to the risk allele
allele.div.2 <- tapply(1:nrow(allele.div), allele.div$pheno, function(x) {
  tt <- allele.div[x, ]
  min.ind <- which.min(tt$beta)
  max.ind <- which.max(tt$beta)
  data.frame(grant=tt$grantham[max.ind]/tt$grantham[min.ind],
            beta=tt$beta[max.ind]/tt$beta[min.ind])
}) %>% do.call(rbind, .) %>% filter_all(all_vars(is.finite(.)))

# add rlm
allele.div.2.mod       <- rlm(beta~grantham, data=data.frame(grantham=allele.div.2$grant, beta=(allele.div.2$beta)), psi=psi.bisquare)
allele.div.2           <- cbind(allele.div.2, predict(allele.div.2.mod, interval='confidence'))
allele.div.2.mod.coeff <- allele.div.2.mod %>% summary %>% .$coefficients %>% .[2, 1]
allele.div.2.mod.p     <- allele.div.2.mod %>% f.robftest(., var="grantham") %>% .$p.value
allele.div.2.mod.cor   <- cor.test(allele.div.2$grant, allele.div.2$beta, method='k')[c('estimate', 'p.value')] %>% unlist
qqplot(allele.div.2.mod$resid, rnorm(10000))

# plot regression for mod allele and risk allele grantham vs beta
p02 <- ggplot(allele.div.2, aes(grant, beta)) +
  geom_point() +
  xlab(expression('grm'['max']*' / '*'grm'['min'])) +
  ylab(expression('beta'['max']*' / '*'beta'['min'])) +
  geom_line(aes(grant, fit), linetype='dashed', size=0.35) +
  geom_label(x=7, y=9, size=3.1, label.size=0, 
             label=paste0("coefficient = ", allele.div.2.mod.coeff %>% signif(., 2), 
                          " (p = ", allele.div.2.mod.p %>% signif(., 2), ")",
                          "\nKendall's tau = ", allele.div.2.mod.cor[1] %>% signif(., 2), 
                          " (p = ", allele.div.2.mod.cor[2] %>% signif(., 2), ")")) +
  ylim(c(2, 10)) +
  #geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3) +
  theme(panel.background=element_blank(), axis.line=element_blank(),
        panel.border=element_rect(fill=NA, size=.2), legend.position=c(.85, .85), legend.title=element_blank(),
        axis.text.x=element_text(angle=0, hjust=1))

cairo_pdf('./results/grantham_reg_17Mar20.pdf', height=4, width=5)
p02
dev.off()



# combined plot
cairo_pdf('./results/genotype_beta_effect_23Mar20.pdf', width=9, height=5.6)
ggarrange(ggarrange(p01, p02, nrow=1, ncol=2, labels='auto'),
          ggtexttable(c(str_wrap("a) Protein-level distances between modulatory alleles. 
                                The plot shows colour-coded density graphs of Grantham's distances for heterozygotic 
                                 genotype associations, haplotype associations and random allele pairs. 
                                 Modulatory alleles in heterozygotic genotypes are significantly different from each other.", 100), 
                        str_wrap("b) Relationship between the risk and modulatory allele distance  
                                and association effect size in heterozygous genotypes. X-axis: Grantham's distance 
                                between risk and modulatory allele (grm) for 
                                maximum / minimum genotype associations. Y-axis: association beta for 
                                maximum / minimum genotype associations.", 100)), 
                      theme=ttheme(tbody.style=tbody_style(hjust=0, x=0, fill=NA, size=10))), 
          nrow=2, ncol=1, heights=c(1, 0.5))
dev.off()


# tbody.style=tbody_style(hjust=0, x=0.1))
ggtexttable(c("a) Protein-level distances between modulatory allele in heterozygotic genotype 
              associations and in haplotype associations.", 
              "b) Relationship between modulatory allele distance to risk allele and association effect 
              size in genotypes. X-axis: Grantham's distance between risk and modulatory allele in heterozygous 
              genotypes for maximum and minimum associations. Y-axis: association beta for maximum and minimum 
              heterozygous genotype associations."), 
            theme=ttheme("blank", tbody.style=tbody_style(hjust=0, x=0.1)))

# plots of relationships of max and min beta and grantham between mod alleles

geno.results.z <- geno.results.z %>% arrange(maxminrel)
geno.results.z$pheno <- factor(geno.results.z$pheno, levels=geno.results.z$pheno) 
cairo_pdf('./results/geno_betarel_17Mar20.pdf', height=7, width=5)
ggplot(geno.results.z, aes( pheno, maxminrel)) +
  geom_col(fill='grey', color=NA) +
  coord_flip() +
  ylab(expression('beta'['max']*' / '*'beta'['min'])) +
  scale_y_continuous(expand=c(0.02, 0)) +
  theme(panel.background=element_blank(), axis.line=element_blank(),
        panel.border=element_rect(fill=NA, size=.2), legend.position=c(.85, .85), legend.title=element_blank(),
        axis.text.x=element_text(angle=0, hjust=1), axis.text.y=element_text(size=5))
dev.off()

geno.results.z <- geno.results.z %>% arrange(grantham)
geno.results.z$pheno <- factor(geno.results.z$pheno, levels=geno.results.z$pheno) 
cairo_pdf('./results/geno_granthamrel_17Mar20.pdf', height=7, width=5)
ggplot(geno.results.z, aes( pheno, grantham)) +
  geom_col(fill='grey', color=NA) +
  coord_flip() +
  ylab("Grantham's distance between mod. alleles") +
  scale_y_continuous(expand=c(0.02, 0)) +
  theme(panel.background=element_blank(), axis.line=element_blank(),
        panel.border=element_rect(fill=NA, size=.2), legend.position=c(.85, .85), legend.title=element_blank(),
        axis.text.x=element_text(angle=0, hjust=1), axis.text.y=element_text(size=5))
dev.off()






##### tests #####


# circular

# https://www.datacamp.com/community/tutorials/sentiment-analysis-R
tmp <- fread('./tmp/prince_dat.tsv', data.table=F) %>%
  filter(decade != "NA") %>% count(decade, charted)
#chorddiag(tmp)
circos.clear()
chordDiagram(tmp)

# circlize
scale100 <-function(x) {
  x <- x + abs(min(x))
  85-(100*(x/max(x)))
}
tmp <- filter(geno.results.anno, pheno=='T1D')[, c('risk.allele', 'mod.allele', 'beta')]
tmp <- arrange(tmp, beta)
tmp$risk.allele <- paste0('risk ', tmp$risk.allele)
tmp.col <- sapply(1-scale100(tmp$beta)/100, function(x) scales::alpha('grey70', x))
tmp.scale <- 1-scale100(tmp$beta)/100
tmp.col.2 <- ifelse(tmp$beta<0, 'forestgreen', 'red')
tmp.col.2 <- c(sapply(seq_along(tmp.scale), function(x) scales::alpha(tmp.col.2[x], tmp.scale[x])), 'black')
names(tmp.col) <- tmp$mod.allele
names(tmp.col.2) <- c(tmp$mod.allele, unique(tmp$risk.allele))
union(tmp[[1]], tmp[[2]])
circos.clear()
circos.par(gap.after=c(10, rep(3, length(unique(tmp$mod.allele))-1), 10))
chordDiagram(tmp, grid.border=NULL, link.border=0, annotationTrack=c("name", "grid"), scale=T,
             grid.col=tmp.col.2, col=tmp.col)
legend('topright', legend='gdds')


drawCircG <- function(dat) {
  tmp <- dat[, c('risk.allele', 'mod.allele', 'beta')]
  tmp <- arrange(tmp, beta)
  tmp$risk.allele <- paste0(dat$locus, '*', tmp$risk.allele)
  tmp.col <- sapply(1-scale100(tmp$beta)/100, function(x) scales::alpha('grey70', x))
  tmp.scale <- 1-scale100(tmp$beta)/100
  tmp.col.2 <- ifelse(tmp$beta<0, 'forestgreen', 'red')
  tmp.col.2 <- c(sapply(seq_along(tmp.scale), function(x) scales::alpha(tmp.col.2[x], tmp.scale[x])), 'black')
  names(tmp.col) <- tmp$mod.allele
  names(tmp.col.2) <- c(tmp$mod.allele, unique(tmp$risk.allele))
  circos.clear()
  circos.par(gap.after=c(10, rep(3, length(unique(tmp$mod.allele))-1), 10), start.degree=50)
  chordDiagram(tmp, grid.border=NULL, link.border=0, annotationTrack=c("name", "grid"), scale=T, preAllocateTracks=0,
               grid.col=tmp.col.2, col=tmp.col, annotationTrackHeight=convert_height(c(1, 2), "mm"))
  title(dat$pheno[1])
}
drawCircG(filter(geno.results.anno, pheno=='T1D'))


cairo_pdf('./results/geno_beta_circle_06Apr20.pdf', height=8, width=8)
op <- par(mfrow=c(3, 3))
sapply(c('T1D', 'M13_RHEUMA', 'RHEUMA_SEROPOS_WIDE', 'M13_ANKYLOSPON', 'L12_PSORIASIS', 'K11_COELIAC', 
  'G6_MS', 'CD2_FOLLICULAR_LYMPHOMA', 'Z21_KIDNEY_TRANSPLANT_STATUS'), function(x) {
    drawCircG(filter(geno.results.anno, pheno==x))
  })
par(op)
dev.off()



drawCircG(filter(geno.results.anno, pheno=='T1D'))
drawCircG(filter(geno.results.anno, pheno=='RHEUMA_SEROPOS_WIDE'))
drawCircG(filter(geno.results.anno, pheno=='CD2_FOLLICULAR_LYMPHOMA'))


# alluvial
library(ggalluvial)
tmp <- filter(geno.results.anno, pheno=='T1D')[, c('risk.allele', 'mod.allele', 'beta')]
tmp$risk.allele[1:3] <- paste0('risk ', tmp$risk.allele[1:3])
tmp <- mutate(tmp, freq=abs(beta/sum(abs(beta))))
ggplot(as.data.frame(tmp), aes(y=freq, axis1=mod.allele, axis2=risk.allele)) +
  geom_alluvium(aes(fill='grey50'), width=1/10)


# ggraph
library(ggraph)
library(igraph)
tmp <- filter(geno.results.anno, pheno=='T1D')[, c('risk.allele', 'mod.allele', 'beta')]
tmp$risk.allele <- paste0('risk ', tmp$risk.allele)
tmp <- tmp[, c(2,1,3)]
colnames(tmp)[1:2] <- c('from', 'to') 
tmp <- igraph::graph_from_data_frame(tmp)
ggraph(tmp, circular = TRUE)

#
  
tmp <- filter(geno.results.anno, pheno=='T1D')[, c('risk.allele', 'mod.allele', 'mod.allele', 'beta')] 
tmp <- spread(tmp, mod.allele.1, beta) %>% data.frame
#tmp[is.na(tmp)] <- 0
tmp <- tmp[, -c(1:2)] %>% as.matrix
dimnames(tmp) <- list(allele=filter(geno.results.anno, pheno=='T1D')$mod.allele, 
                      beta=paste0(filter(geno.results.anno, pheno=='T1D')$risk.allele, ' risk'))

#circos.par(gap.after=10)
circos.par(gap.after=c(rep(10, ncol(tmp)), 20))
chordDiagram(tmp, keep.diagonal=T, scale=T, directional=0,
             row.col=paste0('grey', filter(geno.results.anno, pheno=='T1D')$beta %>% scale100 %>% round))



#
tmp <- lapply(filter(geno.results.anno, pheno=='T1D')$beta, function(x) x*100 %>% round) 
names(tmp) <- filter(geno.results.anno, pheno=='T1D')$mod.allele
BioCircos(genome=tmp, genomeFillColor=c('tomato', 'red'))



# DRB1 on rheuma
HLAdivergence(c('DRB1*04:01', 'DRB1*04:08'), divergences, hla.aligned[['DRB1']], 'grantham')
HLAdivergence(c('DRB1*13:01', 'DRB1*04:08'), divergences, hla.aligned[['DRB1']], 'grantham')

write.table(filter(geno.results.anno, pheno=='T1D') %>% data.frame, './results/T1D_HLAgenotype_all.tsv', sep='\t', row.names=F)

filter(geno.results.anno, pheno=='T1D') %>% data.frame
filter(geno.results.anno, pheno=='T1D_STRICT') %>% data.frame


geno.results.anno[grepl('rheuma', geno.results.anno$pheno, ignore.case=T), 'pheno'] %>% table
geno.results.anno[grepl('rheuma_seropos', geno.results.anno$pheno, ignore.case=T), ]
geno.results.anno$pheno %>% table

cairo_pdf('./results/allele_beta_dist_01Apr20.pdf', height=7.5, width=8.5)
ggarrange(
  plotlist=map(c('T1D', 'M13_RHEUMA', 'RHEUMA_SEROPOS_WIDE', 'M13_ANKYLOSPON', 'L12_PSORIASIS', 'K11_COELIAC', 
                 'G6_MS', 'CD2_FOLLICULAR_LYMPHOMA', 'Z21_KIDNEY_TRANSPLANT_STATUS'), function(x) {
                   tt <- geno.results.anno[geno.results.anno$pheno==x, ] %>% arrange(., beta)
                   tt$beta[tt$beta<0] <- 0.05
                   
                   riskal <- tt$risk.allele[1]#paste(tt$locus[1], tt$risk.allele[1], sep='*')
                   riskcol <- rep('grey60', nrow(tt))
                   riskcol[tt$mod.allele==riskal] <- 'red'
                   tt <- mutate(tt, barcol=riskcol)
                   tt$mod.allele <- factor(tt$mod.allele, levels=unique(tt$mod.allele))
                   ggplot(tt, aes(mod.allele, beta)) +
                     geom_bar(stat='identity', fill=tt$barcol, alpha=0.6) +
                     geom_errorbar(aes(ymin=beta-SEbeta, ymax=beta+SEbeta), width=.2, position=position_dodge(.9), color=riskcol) + 
                     xlab(tt$locus[1]) + ylab(ifelse(x=='M13_ANKYLOSPON', 'beta\n', '')) +
                     ylim(c(-0.3, 3.7)) +
                     facet_wrap(~pheno, scales='free_x') +
                     theme(panel.background=element_blank(), axis.ticks.x=element_blank(),
                           axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8, margin=margin(-6,0,0,0)))
                 }), align='hv'
)
dev.off()


cairo_pdf('./results/allele_beta_dist_flip_02Apr20.pdf', height=8, width=8)
ggarrange(
  plotlist=map(c('T1D', 'M13_RHEUMA', 'RHEUMA_SEROPOS_WIDE', 'M13_ANKYLOSPON', 'L12_PSORIASIS', 'K11_COELIAC', 
                 'G6_MS', 'CD2_FOLLICULAR_LYMPHOMA', 'Z21_KIDNEY_TRANSPLANT_STATUS'), function(x) {
                   tt <- geno.results.anno[geno.results.anno$pheno==x, ] %>% arrange(., beta)
                   #tt$beta[tt$beta<0] <- 0.01
                   tt$SEbeta[tt$beta<0] <- tt$SEbeta[tt$beta<0]*(-1)
                   riskal <- tt$risk.allele[1]#paste(tt$locus[1], tt$risk.allele[1], sep='*')
                   riskcol <- rep('grey75', nrow(tt))
                   riskcol[tt$mod.allele==riskal] <- scales::alpha('red', 0.5)
                   tt <- mutate(tt, barcol=riskcol)
                   riskcol <- rep('grey30', nrow(tt))
                   riskcol[tt$mod.allele==riskal] <- 'red' 
                   tt$mod.allele <- factor(tt$mod.allele, levels=unique(tt$mod.allele))
                   ggplot(tt, aes(mod.allele, beta)) +
                     geom_bar(stat='identity', fill=tt$barcol) +
                     geom_errorbar(aes(ymin=beta, ymax=beta+SEbeta), size=.2, width=0.4, position=position_dodge(.9), 
                                   color=tt$barcol) + 
                     xlab(tt$locus[1]) + ylab(ifelse(x=='CD2_FOLLICULAR_LYMPHOMA', '\nbeta', '')) +
                     #ylim(c(0, 3.7)) +
                     scale_y_continuous(limits=c(ifelse(min(tt$beta+tt$SEbeta)>0, 0, min(tt$beta+tt$SEbeta)), 3.7), expand=c(0, 0)) +
                     coord_flip() +
                     facet_wrap(~pheno, scales='free_x') +
                     theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
                           axis.text.x=element_text(size=9), 
                           axis.text.y=element_text(size=8.5, color=riskcol), axis.line.x=element_line(size=.2),
                           strip.background=element_rect(fill=NA))
                 }), align='hv'
) 
dev.off()




map(c('T1D', 'M13_RHEUMA', 'M13_ANKYLOSPON', 'L12_PSORIASIS', 'K11_COELIAC', 'G6_MS'), function(x) {
  geno.results.anno[geno.results.anno$pheno==x, ] %>% arrange(., beta) })

filter(geno.results, pheno=='M13_RHEUMA')

tmp <- map(c('T1D', 'M13_RHEUMA', 'M13_ANKYLOSPON', 'L12_PSORIASIS', 'K11_COELIAC', 'G6_MS'), function(x) {
  tt <- geno.results.anno[geno.results.anno$pheno==x, ]
  arrange(tt, beta)
}) %>% do.call(rbind, .)
tmp <- unite(tmp, gene, locus, mod.allele, remove=F, sep='-')
tmp$gene <- factor(tmp$gene, levels=unique(tmp$gene))

mapvalues

cairo_pdf('./results/allele_beta_dist_31Mar20.pdf', height=7, width=7)
ggplot(tmp, aes(gene, beta)) +
  geom_bar(stat='identity') +
  facet_wrap(~ pheno, scales='free_x') +
  #scale_x_discrete(labels=(str_split_fixed(tmp$gene, '-', 3)[, 2:3] %>% data.frame %>% 
  #                           unite(., tt, X1, X2, sep='*') %>% .$tt)) +
  theme(panel.background=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(-4,0,0,0)), 
        axis.ticks.x=element_blank())
dev.off()


plot((tmp2$grant), (tmp2$beta))
rlm(beta~grantham, data=data.frame(grantham=tmp2$grant, beta=(tmp2$beta))) %>% summary
qqplot(tmp2.mod$resid, rnorm(10000))
abline(rlm(beta~grantham, data=data.frame(grantham=tmp2$grant, beta=(tmp2$beta))))

plot(tmp2$grant %>% log, tmp2$beta)
plot(tmp2$grant, tmp2$beta %>% log)
plot(tmp2$grant %>% log, tmp2$beta %>% log)
rlm(beta~grantham, data=data.frame(grantham=tmp2$grant, beta=log(tmp2$beta))) %>% summary
glm(beta~grantham, data=data.frame(grantham=log(tmp2$grant), beta=log(tmp2$beta)) %>% filter_all(all_vars(is.finite(.)))) %>% summary

tmp <- tapply(1:nrow(tmp), tmp$pheno, function(x) {
  tt <- tmp[x, ]
  tt$beta <- tt$beta/min(tt$beta)
  #tt$grantham <- tt$grantham/mean(tt$grantham, na.rm=T)
  tt[, c('pheno', 'beta', 'grantham')]
}) %>% do.call(rbind, .)


plot(tmp$grantham, tmp$beta, col=colors()[tmp$pheno %>% factor])
plot(tmp$grantham %>% log, tmp$beta %>% log, col=colors()[tmp$pheno %>% factor])
plot(tmp$grantham, tmp$beta %>% log, col=colors()[tmp$pheno %>% factor])
plot(tmp$grantham %>% log, tmp$beta, col=colors()[tmp$pheno %>% factor])

filter(geno.results.anno, pheno=='T1D') %>% data.frame %>% qplot(grantham, beta, data=.)
filter(tmp, pheno=='H7_IRIDOCYCLITIS') %>% data.frame %>% qplot(grantham, beta, data=., xlim=c(5, 10))
filter(tmp, pheno=='THYROIDITIS') %>% data.frame %>% qplot(grantham, beta, data=.)
filter(tmp, pheno=='ASTHMA_CVDMETABOCOMORB') %>% data.frame %>% qplot(grantham, beta, data=.)
filter(tmp, pheno=='M13_ANKYLOSPON_ICD10') %>% data.frame %>% qplot(grantham, beta, data=., xlim=c(5, 10))


glm(beta~grantham, data=data.frame(grantham=tmp$grantham, beta=log(tmp$beta))) %>% summary
glm(beta~grantham, data=data.frame(grantham=tmp$grantham, beta=(tmp$beta))) %>% summary
glm(beta~grantham, data=data.frame(grantham=(tmp$grantham), beta=(tmp$beta))) %>% summary


#Z.test.2(geno.results[2, 'beta'], geno.results[3, 'beta'], geno.results[2, 'SEbeta'], geno.results[3, 'SEbeta'])
geno.results.z <- tapply(1:nrow(geno.results.anno), geno.results.anno$pheno, function(x) {
  tmp <- geno.results.anno[x, ]
  tmp <- mutate(tmp, secoef=abs(beta)/SEbeta)
  #tmp <- filter(tmp, secoef>1)
  ind.max <- which.max(tmp$beta)
  ind.min <- which.min(tmp$beta)
  out <- Z.test.2(tmp[ind.max, 'beta'], tmp[ind.min, 'beta'], tmp[ind.max, 'SEbeta'], tmp[ind.min, 'SEbeta'])
  data.frame(pheno=tmp$pheno[1], locus=tmp$locus[1], risk.allele=tmp$risk.allele[1],
             mod.allele.max=tmp$mod.allele[ind.max], beta.max=tmp$beta[ind.max], beta.max.SE=tmp$SEbeta[ind.max],
             mod.allele.min=tmp$mod.allele[ind.min], beta.min=tmp$beta[ind.min], beta.min.SE=tmp$SEbeta[ind.min],
             p.value=out)
}) %>% do.call(rbind, .)

