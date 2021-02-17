library(data.table)
library(ggplot2)
library(scales)
## work in the tmp directory
## save figures to exportDir
exportDir = ''
### MIDPOINT IS SET TO 0.333 BECAUSE WE HAVE 3 BODYSITES SO 0.33 PROB PER SITE IF RANDOMLY DISTRIBUTED 
EC_asv = 'ASV_1:Escherichia/Shigella'
SC_asv = "ASV_3:Staphylococcus"
LC_asv = "ASV_2:Lactobacillus"



df = fread('ClustersInBodySites.csv')
df = df[, SampleType:= factor(SampleType, levels=c('Gut','Oral cavity','Skin'))]
df <- df[order(SampleType),]
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == EC_asv,'EC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == SC_asv,'SC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == LC_asv,'LC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == 'Intermediate','IC')
ggplot(df, aes(x=SampleType, y=ClusterASV, fill=value)) + geom_tile(size=.1) + 
  scale_fill_gradient2('Samples in cluster (%)',labels=percent, midpoint = .33) +
  theme(aspect.ratio = 1.4, axis.title = element_text(size=16), axis.text = element_text(size=16), 
        legend.text= element_text(size=16), legend.title = element_text(size=14))+
  labs(y='Cluster', x='Body site') 

ggsave(paste0(exportDir,'ClusterInBodysites_heatmap.svg'), height = 6, width = 6)
