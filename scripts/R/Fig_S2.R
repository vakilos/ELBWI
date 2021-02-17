library(ggplot2)
library(data.table)
## work in the tmp directory
## save figures to exportDir
exportDir = ''


df = fread('masterTable.tsv')
df$SampleType = replace(df$SampleType, df$SampleType == 'Oral_cavity','Oral cavity')
brown = '#653700'
pink = '#ff796c'
orange = '#f1c50d'

ggplot(df, aes(PC1, PC2,color=SampleType)) + geom_point(size=5, alpha=0.6) +
  scale_color_manual("Body site",values=c('Gut'=brown, 'Oral cavity'=pink, 'Skin'=orange)) + 
  xlab("PC1 (28%)") + ylab('PC2 (22%)') + theme_bw() + theme(axis.title= element_text(size =16)) +
  scale_shape('Timepoint')+ guides(color=guide_legend(order=1)) + ggtitle("PCoA on Bray-Curtis Dissimilarities") +
  theme(legend.key = element_blank(), legend.title=element_text(size=12), 
        legend.text=element_text(size=10))




ggsave(paste0(exportDir,'PcoA_bodysite_R.svg'), width=10, height=8)
ggsave(paste0(exportDir,'PcoA_bodysite_R.png'), dpi=200)






### Run suppl_figures.py first 
df = fread('WithinPrev_Prev_ToR.csv')

df
### 
ggplot(df, aes(x=Prevalence,y=InPrev, size=MedianAbundance)) + geom_point(alpha=.3) + 
  geom_text(data=subset(df, InPrev > 0.6 & Prevalence < 0.12 & MedianAbundance >= .05), aes(Prevalence,InPrev,label=Taxon), size=3,vjust=-1, hjust=-0.1) +
  geom_point(data=subset(df, InPrev >0.6 & Prevalence < 0.12 & MedianAbundance >= .05),aes(Prevalence,InPrev, color='red')) + 
  theme(aspect.ratio = 0.9,axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  labs(y='Max Prevalence within Patient (%)' ,x ='Dataset Prevalence (%)') + 
  scale_color_discrete('Median relative abundance > 10%\nPrevalence < 10%\nMax prevalence within patient > 50%') + scale_size('Median relative abundance')

ggsave(paste0(exportDir,'Fig_S2.svg'), width = 7, height=5)

    
