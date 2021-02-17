
library(data.table)
library(ggplot2)
library(gridExtra)
library(svglite)
library(ggalt)
library(MASS)
library(ggthemes)
library(ggpubr)
library(grid)

## work in the tmp directory
## save figures to exportDir
exportDir = ''



#### Copy ASV IDs from the output of preterm_paper_analysis.py
EC_asv = 'ASV_1:Escherichia/Shigella' ### orange color for cluster
SC_asv = "ASV_3:Staphylococcus" ### blue color for cluster
LC_asv = "ASV_2:Lactobacillus" ### green color for cluster

ma = fread('masterTable.tsv')

ma$ClusterASV = replace(ma$ClusterASV, ma$ClusterASV == EC_asv,'EC')
ma$ClusterASV = replace(ma$ClusterASV, ma$ClusterASV == SC_asv,'SC')
ma$ClusterASV = replace(ma$ClusterASV, ma$ClusterASV == LC_asv,'LC')
ma$ClusterASV = replace(ma$ClusterASV, ma$ClusterASV == 'Intermediate','IC')


ggplot(data = ma, mapping = aes(x = PC1, y = PC2)) +
  geom_pointdensity(adjust=4) +
  scale_color_viridis()
p0 = kde2d(ma$PC1, ma$PC2, n=100)
pdf(paste0(exportDir,'fig3_kde2d.pdf'), height = 10, width = 10)
filled.contour(p0, xlab = "PC1",ylab = "PC2") 
dev.off()
p1 = ggplot(ma, aes(x=PC1,y=PC2,color=ClusterASV)) + 
  geom_encircle(inherit.aes = TRUE, fill='#f3f4f5', alpha=.3, show.legend = FALSE) +
  geom_point(aes(size=ClusterASV), alpha=0.8)+
  scale_size_manual(values=c(5, 5,5,5))+ theme_bw()+theme(legend.position = 'none', axis.title  = element_text(size=16)) + 
  labs(x='PC1 [28 %]', y='PC2 [22 %]')+
  scale_color_manual(values=c("EC"='#fca481',"LC"='#7dc7b0',
                              "SC"= '#919ebb', 'IC'='#eba1cf')) +
  scale_fill_manual(values=c("EC"='#fca481',"LC"='#919ebb',
                             "SC"= '#7dc7b0', 'IC'='#eba1cf')) +
  annotate('text',label="SC",x=-0.25,y=-0.15, size=6) + annotate('text',label="EC",x=0.27,y=0.21, size=6) +
  annotate('text',label="LC",x=0.35,y=-0.25, size=6) + annotate('text',label="IC",x=-0.03,y=0.2, size=6)
p1
ggsave(paste0(exportDir,'PCOA_clusters.svg'), height=8, width = 9)

df = fread('DominantASVperCluster.tsv')
df
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == EC_asv,'EC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == SC_asv,'SC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == LC_asv,'LC')
df$ClusterASV = replace(df$ClusterASV, df$ClusterASV == 'Intermediate','IC')

df
df = df[, ClusterASV:= factor(ClusterASV, levels=c("EC","LC","SC", 'IC'))]
levels(df$ClusterASV)


ec = ggplot(subset(df, df$variable == EC_asv), aes(x=ClusterASV, y=value, fill=ClusterASV)) + geom_boxplot() +
  scale_fill_manual(name='Cluster' , labels=c('EC','LC','SC','IC') , 
                    values=c("EC"='#fca481',"LC"='#7dc7b0',
                             "SC"= '#919ebb', 'IC'='#eba1cf'))+ labs(x='', y='Relative Abundance',title=EC_asv)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=18), plot.title = element_text(size=18,hjust=0.5),
         aspect.ratio = 1, legend.position = 'none')+ scale_y_continuous(breaks=c(0,0.5,1.0))+
  stat_compare_means(comparisons = list(c('EC','LC'), c('EC','IC'), c('EC','SC')), test = 'wilcox.test',label = "p.signif", size=6, vjust=0.5)
ec
ggsave(paste0(exportDir,'Fig2_EC.svg'), width=5, height=5 )
lac = ggplot(subset(df, df$variable == LC_asv), aes(x=ClusterASV, y=value, fill=ClusterASV)) + geom_boxplot() +
  scale_fill_manual(name='Cluster' , labels=c('EC','LC','SC','IC') , 
                    values=c("EC"='#fca481',"LC"='#7dc7b0',
                             "SC"= '#919ebb', 'IC'='#eba1cf'))+ labs(x="Cluster",y="", title = LC_asv)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        plot.title = element_text(size=18,hjust = 0.5), 
        aspect.ratio = 1, legend.position = 'none')+ scale_y_continuous(breaks=c(0,0.5,1.0))+
  stat_compare_means(comparisons =list(c('EC','LC'), c('LC','IC'), c('LC','SC')) , test = 'wilcox.test',label = "p.signif", size=6, vjust=0.5)
lac
ggsave(paste0(exportDir,'Fig2_LC.svg'), width=5, height=5 )
sta = ggplot(subset(df, df$variable == SC_asv), aes(x=ClusterASV, y=value, fill=ClusterASV)) + geom_boxplot() +
  scale_fill_manual(name='Cluster' , labels=c('EC','LC','SC','IC') , 
                    values=c("EC"='#fca481',"LC"='#7dc7b0',
                             "SC"= '#919ebb', 'IC'='#eba1cf'))+ labs(x='',y='',title = SC_asv)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        plot.title = element_text(size=18,hjust=0.5),
        aspect.ratio = 1, legend.position = 'none')+ scale_y_continuous(breaks=c(0,0.5,1.0))+
  stat_compare_means(comparisons = list(c('SC','LC'), c('SC','IC'), c('EC','SC')), test = 'wilcox.test',label = "p.signif", size=6, vjust=0.5)
ggsave(paste0(exportDir,'Fig2_SC.svg'), width=5, height=5 )

grid.newpage()
g = grid.arrange(ec, lac, sta, nrow=1)

svglite(paste0(exportDir,'ASVClusterAbundance.svg'), width = 12, height = 6)
grid.draw(g)
dev.off()


plot = ggplot(subset(df, df$variable == SC_asv), aes(x=ClusterASV, y=value, fill=ClusterASV)) + geom_boxplot() +
  scale_fill_manual(name='Cluster' , labels=c('EC','LC','SC','IC') , 
                    values=c("EC"='#fca481',"LC"='#7dc7b0',
                             "SC"= '#919ebb', 'IC'='#eba1cf'))+
  theme(axis.text = element_text(size=12), axis.text.x = element_blank(),
        axis.title = element_text(size=18),
        strip.text = element_text(size=14), aspect.ratio = 0.9)
legend <- cowplot::get_legend(plot)
grid.newpage()
svglite(paste0(exportDir, 'legend_ASVclusterAbundance.svg'))
grid.draw(legend)
dev.off()





