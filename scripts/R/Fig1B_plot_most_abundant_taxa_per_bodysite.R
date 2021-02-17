library(ggplot2)
#library(ggridges)
library(data.table)
library(dplyr)
library(viridis)
library(svglite)

## work in the tmp directory
## save figures to exportDir
exportDir = ""

df = fread('mostabuperbodysitephylum.tsv')

levels = df[, .(Relative_abundance = mean(Relative_abundance)), by = .(Taxonomy)][order(Relative_abundance), Taxonomy]
df[ ,Taxonomy := factor(Taxonomy, levels)]

phylum = ggplot(df, aes(y=Relative_abundance, x= Taxonomy)) + geom_boxplot(fill='grey') + 
  coord_flip() + facet_wrap(~SampleType) + 
  theme(strip.text.x = element_text(size = 12), aspect.ratio = 1, legend.position = 'none',panel.spacing = unit(2, "lines"), 
        axis.text.x = element_text(angle=0,size=10, vjust=0.4), axis.text.y=element_text(size=12), axis.title = element_text(size=14))+
  labs(y='Relative abundance', x='Phylum') 

phylum
ggsave(file=paste0(exportDir,'MostAbundantPhylum.svg'), plot=phylum, height=9, width=7, dpi=200)





