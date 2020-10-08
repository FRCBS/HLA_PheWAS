
## libs
sapply(c('tidyverse','data.table','ggpubr','viridis','rasterpdf','mutoss', 'Biostrings', 'MASS', 
         'sfsmisc', 'circlize', 'BioCircos', 'ghibli', 'RColorBrewer', 'SPAtest', 'ggstance', 'ggrepel'),
       library, character.only=T)

options(stringsAsFactors=F)


## --------------------------------------------------------
## functions
## --------------------------------------------------------

# beta difference Z test (two-way)
Z.test.2 <- function(m1, m2, se1, se2) { # means and sd's
  zstat <- (m1-m2) / sqrt( (se1^2) + (se2^2) )
  pnorm(abs(zstat), lower.tail=F)*2
}


# HLA notation formatting
formHLA <- function(x) {
  out <- unite(str_split_fixed(x, '\\.', 2) %>% data.frame, hla, X1, X2, sep='*')$hla
  out <- unite(str_split_fixed(out, '\\.', 2) %>% data.frame, hla, X1, X2, sep=':')$hla
  gsub(':$', '', out)
}

# data read-in and MHC manhattan plotting for SNP summary stats
multiManhattan <- function(phe.ind, ann='') {

  # read in a data file and reformat
  dat <- lapply(ls.snps[phe.ind], function(i) {
    dat     <- fread(i, data.table=F)
    colnames(dat) <- c('chrom',	'pos', 'ref', 'alt', 'rsids', 'nearest_genes', 'pval',	'beta',	'sebeta', 'maf', 'maf_cases',	
                       'maf_controls', 'n_hom_cases', 'n_het_cases', 'n_hom_controls', 'n_het_controls')
    dat$pos <- dat$pos/1e6 
    dat     <- mutate(dat, logP=log10(pval)*(-1))
    tmp.nm  <- gsub('./data/sumstats_MHC_finngen_R5_all/MHC_finngen_R5_', '', i, fixed=T)
    dat     <- data.frame(dat, pheno=tmp.nm)
  }) %>% do.call(rbind, .)
 
  # MHC region gene names and colors for plotting 
  dat.col.vec <- rep('no', nrow(dat))
  sapply(1:nrow(mhc.genes.selected), function(i) {
    dat.col.vec[
      which(dat$pos >= mhc.genes.selected$start_position[i] & dat$pos <= mhc.genes.selected$end_position[i])
      ] <<- 'gene' 
  }) %>% invisible
  dat <- data.frame(dat, gene=dat.col.vec %>% factor(levels=dat.col.vec %>% unique))
  
  # annotation data frame for top SNP
  dat.snp <- tapply(1:nrow(dat), dat$pheno, function(x) {
    tmp <- filter(dat[x, ], pval==min(pval))
    tmp$nearest_genes <- paste0('(', tmp$nearest_genes, ')')
    tmp
  }) %>% do.call(rbind, .)
  dat.snp <- unite(dat.snp, anno, rsids, nearest_genes, sep=' ')
  
  # draw and return a ggplot
  ggplot(dat, aes(pos, logP)) +
    geom_point(size=0.99, fill='#a0c4e2', color='#5aadf2', shape=21, alpha=0.7, stroke=0.15) +
    geom_point(aes(pos, logP), data=filter(dat, gene=='gene'), fill='royalblue1', color='royalblue3', 
               size=0.99, shape=21, stroke=0.15) +
    geom_hline(yintercept=(log10(5e-8)*(-1)), linetype='dashed', size=0.2, alpha=0.9, color='black') +
    ylab(expression('-log'[10]*'('*italic(p)*')')) +
    scale_x_continuous(limits=c(29.7, 33.12), breaks=mhc.genes.selected$gene_mean, 
                       labels=mhc.genes.selected$wikigene_name, 
                       position='bottom', name='gene',  minor_breaks=NULL, guide=guide_axis(check.overlap=F),
                       sec.axis=sec_axis(~., name='chr 6 position (Mb)'),
                       expand=expansion(0.01)) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
    facet_wrap(~pheno, ncol=1, strip.position='right', scale='free_y') +
    geom_label_repel(aes(pos, logP, label=anno), data=dat.snp, segment.size=0.2, nudge_x=0.3, nudge_y=-0.3, 
                     size=2.5, min.segment.length=0) +
    geom_text(data=ann, mapping=aes(x=29.7, y=Inf, label=label), hjust=0, vjust=2, size=3.2) +
    theme(panel.background=element_blank(), panel.border=element_rect(size=0.2, fill=NA), legend.position='none',
          strip.background=element_rect(fill='white'), axis.text.x.bottom=element_text(angle=45, hjust=1, size=7),
          axis.text.x.top=element_text(size=8), strip.text=element_text(size=7),
          axis.line.x=element_line(size=.3), axis.line.y=element_line(size=.3)) 
}


