library(ggplot2)
library(data.table)
library(gtable)
library(ggpubr)
library(grid)
library(svglite)
library(forcats)
#### work in the tmp directory

## export Figures to this directory
exportDir = ''


df=fread('prev_toR.tsv')
we = fread('mostabuperbodysitegenus.tsv')
### Based on the output of preterm_paper_analysis
we[we== "ASV_29"]<-"unclassified_Enterobacteriaceae"
df[df== "ASV_29"]<-"unclassified_Enterobacteriaceae"

showTaxa=c('Staphylococcus', 'Lactobacillus', 
'Escherichia/Shigella', 'Enterococcus', 'Proteus', 
'Corynebacterium', 'Ureaplasma', 'Bacteroides', 'Prevotella', 
'Fusobacterium', 'Streptococcus', 'unclassified_Enterobacteriaceae', 'Enterobacter', 'Lachnoanaerobaculum')


#### ANOVA
staph_anova = aov(Relative_abundance ~ SampleType, data=subset(we, we$Taxonomy == 'Staphylococcus'))
esche_anova = aov(Relative_abundance ~ SampleType, data=subset(we, we$Taxonomy == 'Escherichia/Shigella'))
lact_anova = aov(Relative_abundance ~ SampleType, data=subset(we, we$Taxonomy == 'Lactobacillus'))
summary(staph_anova)
summary(esche_anova)
summary(lact_anova)
TukeyHSD(staph_anova)
TukeyHSD(lact_anova)
TukeyHSD(esche_anova)
unique(we$Taxonomy)

we = subset(we, Taxonomy %in% showTaxa)
df = subset(df, Taxon %in% showTaxa)
we = we[, Taxonomy:= factor(Taxonomy, rev(showTaxa))]
df = df[, Taxon := factor(Taxon, rev(showTaxa))]
we


for (s in c('Gut', "Oral cavity",'Skin')){
  name = toString(s)
  p2 = ggplot(subset(we,we$SampleType == s), aes(y=Relative_abundance, x= Taxonomy)) + geom_boxplot(fill='grey') + 
    coord_flip() + 
    theme(strip.text.x = element_text(size = 12), aspect.ratio = 2, legend.position = 'none',panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(angle=0,size=10, vjust=0.4), axis.text.y=element_text(size=12), axis.title = element_text(size=14))+
    labs(title = name, y='Relative abundance', x='Genus') 
  
  p3 = ggballoonplot(subset(df, df$Bodysite == s), fill='value') + scale_fill_viridis_c(option = "C")+ 
    coord_flip() + theme(legend.position = 'none', axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid.major = element_blank())
  g2 <- ggplotGrob(p2)
  g3 <- ggplotGrob(p3)
  g <- cbind(g2, g3, size='first')
  #g$heights <- unit.pmax(g2$heights, g3$heights)
  grid.newpage()
  svg(paste0(exportDir,'AbundancePrevelance_',name,'.svg'))
  grid.draw(g)
  dev.off()
  }


plot = ggballoonplot(df, fill='value') + scale_fill_viridis_c(name='Prevalence',option = "C")+ coord_flip() +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid.major = element_blank())
legend <- cowplot::get_legend(plot)
grid.newpage()
svglite(paste0(exportDir,'abundance_prevalence_legend.svg'))
grid.draw(legend)
dev.off()
      
  
