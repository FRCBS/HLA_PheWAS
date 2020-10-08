

# libs and functions
source('./src/functions.R')



## --------------------------------------------------------
## plot the PheWAS summary
## --------------------------------------------------------

# top phenotypes
p.hla.phewas <- ggplot(results.anno.2.1, aes(logP2, DESC)) +
  geom_jitter(color=results.anno.2.1$col, fill=results.anno.2.1$col,
              size=logp2size(results.anno.2.1$logP2), alpha=0.5, shape=21, height=0.4) +
  scale_x_continuous(limits=c(0, 2.7), breaks=c(0, 1, 2), expand=c(0, 0)) +
  xlab(expression('log'[10]*'(-log'[10]*'('*italic(p)*')+1)')) + ylab('phenotype category') +
  geom_vline(xintercept=log10((log10(5e-8)*(-1))+1), linetype='dashed', size=0.3, alpha=0.8) +
  geom_vline(xintercept=log10((log10(filter(results.anno, adjusted.p<0.01)$p.value %>% max)*(-1))+1), 
             linetype='dashed', size=0.3, alpha=0.8) +
  facet_wrap(~ gene, nrow=1, ncol=7) +
  theme(panel.background=element_blank(), axis.ticks.y=element_blank(), 
        axis.text.x=element_text(size=7, angle=0, hjust=0.5), axis.text.y=element_text(size=8.6), 
        panel.border=element_rect(fill=NA, size=.3),
        axis.title.y=element_text(color='black'), legend.position='none',
        strip.background=element_rect(fill=NA)) 



## --------------------------------------------------------
## plot comparison with other HLA phewas
## --------------------------------------------------------

p.studies.comp <- ggplot(studies.comp, aes(odds_ratio_study, odds_ratio_y)) +
  geom_point(aes(fill=phenotype, shape=study, group=group, color=phenotype), alpha=0.95, size=4) +
  geom_label_repel(aes(odds_ratio_study, odds_ratio_y, label=snp), segment.size=0.2, nudge_x=0.9, nudge_y=-0.9,
                   data=studies.comp[c(163, 168, 174, 148, 153, 159), ], 
                   size=2.3, segment.alpha=0.3, color='grey35', alpha=0.8) +
  geom_label_repel(aes(odds_ratio_study, odds_ratio_y, label=snp), segment.size=0.2, nudge_x=-0.7, nudge_y=0.7,
                   data=studies.comp[c(116, 118, 85, 2, 125), ], 
                   size=2.3, segment.alpha=0.3, color='grey35', alpha=0.8) +
  scale_shape_manual(values=c(24, 21, 23)) +
  scale_color_manual(values=rep('white', studies.comp$phenotype %>% unique %>% length)) +
  scale_fill_manual(guide=guide_legend(override.aes=list(shape=22)),
    values=c(brewer.pal(8, name='Set2'), brewer.pal(8, name='Dark2'))[1:(studies.comp$phenotype %>% unique %>% length)]) +
  xlab('odds ratio (previous studies)') + ylab('odds ratio (present study)') +
  xlim(c(0, 7.1)) + ylim(c(0, 6.1)) +
  geom_hline(yintercept=1, size=0.2, linetype='dashed', color='grey40') +
  geom_vline(xintercept=1, size=0.2, linetype='dashed', color='grey40') +
  theme(panel.background=element_blank(),  legend.key=element_blank(),
        axis.text.x=element_text(size=8, angle=0, hjust=1), axis.text.y=element_text(size=8.5), 
        axis.line.x=element_line(size=.3), axis.line.y=element_line(size=.3), 
        axis.title.y=element_text(color='black'), legend.position='right',
        strip.background=element_rect(fill=NA)) 



## --------------------------------------------------------
## plot replication correlations
## --------------------------------------------------------

p.rep.corr <- ggplot(rbind(data.frame(filter(results.anno, adjusted.p<0.01, beta<2.5)[c('beta', 'beta.replication')], group='allele'),
                           data.frame(filter(diplo.results.anno, adjusted.p<0.01)[c('beta', 'beta.replication')], group='diplotype'),
                           data.frame(filter(haplo.results.anno, adjusted.p<0.01)[c('beta', 'beta.replication')], group='haplotype')),
                     aes(x=beta, y=beta.replication)) +
  geom_point(shape=21, color='grey10', alpha=0.4, size=0.3) +
  geom_smooth(method='lm', formula=y~x, size=0.3, alpha=0.7, se=F) +
  stat_cor(method="pearson", aes(label=..r.label..), size=3) +
  xlab('beta (discovery cohort, FDR<0.01)') + ylab('beta\n(replication cohort)') +
  facet_wrap(~group, scales='free') +
  theme(strip.background=element_rect(fill=NA, color='white'), strip.text=element_text(size=12),
        panel.border=element_blank(), panel.background=element_blank(), axis.line.x=element_line(size=0.3),
        axis.line.y=element_line(size=0.3), axis.text.y=element_text(size=8.5), axis.text.x=element_text(size=8.5),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())



## --------------------------------------------------------
## plot shared alleles by tag
## --------------------------------------------------------

# plot the number of independent categories per HLA 
p.allele.count <- ggplot(results.cond.f.tags.2, aes(HLA, DESC, fill=factor(hla_n))) +
  geom_tile(color='black', alpha=0.9) +
  scale_fill_manual(values=viridis(results.cond.f.tags.2$hla_n %>% unique %>% length), 
                    name='independent\ncategories') +
  theme_minimal() +
  ylab('phenotype category') + xlab('HLA') +
  theme(legend.position='top', legend.direction="horizontal", 
        legend.key.size=unit(0.75, 'line'), legend.margin=margin(t=0, r=0, b=0, l=0),
        legend.title=element_text(size=10), legend.text=element_text(size=8),
        axis.text.x=element_text(angle=45, hjust=1, size=5.5, margin=margin(-1, 0,0,0)), 
        axis.text.y=element_text(angle=0, hjust=1, size=8.5), 
        axis.title.y=element_text(margin=margin(t=0, r=14, b=0, l=0)),
        panel.border=element_rect(size=0.3, fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white", colour="white"),
        axis.line.x=element_line(size=0), axis.line.y=element_line(size=0)) +
  guides(fill=guide_legend(nrow=1, byrow=T, label.position='top', title.position='left', reverse=T))

# infectious/autoimmune
p.allele.count.2 <- ggplot(results.cond.f2, aes(covpheno, longname, fill=beta)) +
  geom_tile(color='black') +
  xlab('autoimmune disease phenotype') + ylab('infectious disease\nphenotype') +
  facet_wrap(~HLA, nrow=3, drop=F) +
  theme_minimal() +
  scale_fill_gradient2(low="seagreen", high="red", na.value="grey") +
  theme(axis.text.x=element_text(angle=24, hjust=1, size=6), axis.text.y=element_text(size=8), 
        axis.title.y=element_text(margin=margin(t=0, r=38, b=0, l=0)),
        panel.border=element_rect(size=0.3, fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white", colour="white"),
        legend.key.size=unit(0.75, 'line'), legend.margin=margin(t=0, r=0, b=0, l=0),
        legend.title=element_text(size=10), legend.text=element_text(size=8)) 



## --------------------------------------------------------
## Figures
## --------------------------------------------------------

# Figure 1
jpeg('./results/figures/Figure1.jpg', height=4, width=10, units='in', res=600)                           
p.hla.phewas
dev.off()

# Figure 2
jpeg('./results/figures/Figure2.jpg', width=8, height=7, res=600, unit='in')
ggarrange(ggarrange(ggplot()+theme_void(), p.studies.comp, widths=c(0.034, 1)), 
          ggplot()+theme_void(),
          ggarrange(p.rep.corr, ggplot()+theme_void(), widths=c(1, 0.33)), 
          ncol=1, nrow=3, heights=c(1, 0.01, 0.45), labels=c('A', '', 'B'), font.label=list(size=15))
dev.off()

# Figure 3
jpeg('./results/figures/Figure3.jpg', height=6.5, width=10, units='in', res=600)                           
ggarrange(p.allele.count, 
          p.allele.count.2, 
          ncol=1, nrow=2, labels=c('A', 'B'), heights=c(1, 1), widths=1, font.label=list(size=15))
dev.off()

# Figure 4
jpeg('./results/figures/Figure4.jpg', height=10, width=9, unit='in', res=600)
p.diplo.betas
dev.off()

# Figure 5
jpeg('./results/figures/Figure5.jpg', height=6.5, width=9.5, res=600, unit='in')
p.haplo.betas
dev.off()

# Figure 6
jpeg('./results/figures/Figure6.jpg', res=600, unit='in', width=9, height=10)
ggarrange(p.snps01, p.snps02, heights=c(3, 2), labels=c('A', 'B'), ncol=1, nrow=2)
dev.off()


