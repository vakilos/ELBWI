library(gplots)
library(ggpubr)
library(data.table)
## work in the tmp directory
## save figures to exportDir
exportDir = ''

af = fread('MonoFreq_toR.csv')
bf = fread('MonoFreq_bpd_toR.csv')
af
xtabs(~BPD + SampleType, data = bf )
xtabs(~IVH + SampleType, data = af)
ggplot(af, aes(y=MonoFreq, x=IVH, fill=IVH)) + geom_violin(width=.4)+ geom_jitter(size=0.2) + facet_wrap(~SampleType) + 
  scale_fill_manual(values = c('#d39b9b', '#d16666')) + labs(y = 'Frequency of mono-dominated community state', x='IVH diagnosis') + 
  stat_compare_means(method='wilcox.test', label='p.signif',label.x = 1.3 , size= 8, vjust = 0.1,hjust=-1, hide.ns = TRUE ) + 
  theme(strip.text = element_text(size=16), axis.title = element_text(size=16),
        axis.text.x = element_text(size=14), legend.position = 'none', aspect.ratio = 2)

  ggsave(paste0(exportDir,'MonoDominantClusterFrequencyIVH.svg'))
  ## Wilcoxon p-adj = 0.044, p= 0.041, p = 0.023
compare_means(method='wilcox.test', data=subset(af, af$SampleType == 'Gut'), MonoFreq~IVH)
compare_means(method='wilcox.test', data=subset(af, af$SampleType == 'Oral cavity'), MonoFreq~IVH)
compare_means(method='wilcox.test', data=subset(af, af$SampleType == 'Skin'), MonoFreq~IVH)

############################################################################
### for BPD
ggplot(bf, aes(y=MonoFreq, x=BPD, fill=BPD)) + geom_violin(width=.4)+ geom_jitter(size=0.2) + facet_wrap(~SampleType) + 
  scale_fill_manual(values = c('#d39b9b', '#d16666')) + labs(y = 'Frequency of mono-dominated community state', x='BPD diagnosis') + 
  stat_compare_means(method='wilcox.test', label='p.signif',label.x = 1.3 , size= 8, vjust = 0.1,hjust=-1, hide.ns = TRUE ) + 
  theme(strip.text = element_text(size=16), axis.title = element_text(size=16),
        axis.text.x = element_text(size=14), legend.position = 'none', aspect.ratio = 2)

ggsave(paste0(exportDir,'MonoDominantClusterFrequencyIVH.svg'))
## Wilcoxon p-adj = 0.044, p= 0.041, p = 0.023
compare_means(method='wilcox.test', data=subset(bf, af$SampleType == 'Gut'), MonoFreq~BPD)
compare_means(method='wilcox.test', data=subset(bf, af$SampleType == 'Oral cavity'), MonoFreq~BPD)
compare_means(method='wilcox.test', data=subset(bf, af$SampleType == 'Skin'), MonoFreq~BPD)








