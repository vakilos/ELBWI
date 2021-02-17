library(data.table)
library(ggplot2)
library(ggpubr)
## work in the tmp directory
## save figures to exportDir

#### get ASV IDs from preterm_paper_analysis

EC_asv = 'ASV_1:Escherichia/Shigella'
SC_asv = "ASV_3:Staphylococcus"
LC_asv = "ASV_2:Lactobacillus"
exportDir = ''

df = fread('masterTable.tsv')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == EC_asv,'EC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == SC_asv,'SC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == LC_asv,'LC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == 'Intermediate','IC')
df = df[, ClusterASV:= factor(ClusterASV, levels=c("EC","LC","SC", 'IC'))]
df <- df[order(ClusterASV),]

ggplot(df, aes(x=ClusterASV, y=AlphaShannon, fill=ClusterASV)) + geom_boxplot() + labs(x='Cluster', y= 'Shannon index')+ 
  theme(aspect.ratio = 1.7, axis.title = element_text(size=18), axis.text = element_text(size=18)) +
  scale_fill_manual(name='Cluster' , labels=c('EC','LC','SC','IC') , 
                    values=c("EC"='#fca481',"LC"='#7dc7b0',
                             "SC"= '#919ebb', 'IC'='#eba1cf')) +
    stat_compare_means(comparisons=list( c('IC','EC'),c('IC','SC'), c('IC','LC')), size=8,test= 'wilcox.test',label = "p.signif", vjust=0.5)

ggsave(paste0(exportDir, 'Shannon_between_clusters.svg'))

compare_means(AlphaShannon ~ ClusterASV, data=df, method='wilcox.test')

