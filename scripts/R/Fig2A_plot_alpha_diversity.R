library(data.table)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
library(ggExtra)
library(gtable)
library(gridExtra)
library(svglite)
library(ggpubr)
## work in the tmp directory
## save figures to exportDir
exportDir = ''

df = fread('CountTable.tsv', key='SampleID')
mf = fread('MetaTable.tsv', key='SampleID')
mf$SampleType = replace(mf$SampleType, mf$SampleType == 'Oral_cavity','Oral cavity')
mf$SampleType
brown = '#653700'
pink = '#ff796c'
orange = '#f1c50d'


grid_arrange_shared_legend <-
function(myList,
         ncol = length(myList),
         nrow = 1,
         position = c("bottom", "right")) {

  plots <- myList
  position <- match.arg(position)
  g <-
    ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x)
    x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x)
    x + theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(
    position,
    "bottom" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight)
    ),
    "right" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 1,
      widths = unit.c(unit(1, "npc") - lwidth, lwidth)
    )
  )

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}
af = merge(df,mf, by='SampleID')

shan = ggplot(af, aes(Age, AlphaShannon, colour = SampleType)) + geom_point()+
  scale_color_manual(name='Body Site',values=c('Gut'=brown, 'Oral cavity'=pink, 'Skin'=orange))+ ylab('Shannon index')+ xlab('DoL')+
  scale_fill_manual(name='Body Site',values=c('Gut'=brown, 'Oral cavity'=pink, 'Skin'=orange))+
  geom_smooth(method='lm') + facet_wrap(~SampleType)+ xlab('') + stat_cor(method = 'pearson', label.x = 0.1,label.y=2.5) +
  theme(strip.text.x = element_text(size = 14), axis.title= element_text(size =18), axis.text = element_text(size=12),
        legend.position = 'none')+theme(aspect.ratio = .8) 
shan
#f = ggMarginal(plot, groupColour = TRUE, groupFill = TRUE)
chao = ggplot(af, aes(Age, AlphaChao, colour = SampleType)) + geom_point()+
  scale_color_manual(name='Body Site',values=c('Gut'=brown, 'Oral cavity'=pink, 'Skin'=orange))+ ylab('Chao index')+xlab('DoL')+
  scale_fill_manual(name='Body Site',values=c('Gut'=brown, 'Oral cavity'=pink, 'Skin'=orange))+
  geom_smooth(method='lm') + facet_wrap(~SampleType)+ xlab('') +  stat_cor(method = 'pearson', label.x = 0.1,label.y=20) +
  theme(strip.text.x = element_text(size = 14),axis.title= element_text(size =18),axis.text = element_text(size=12),
        legend.position = 'none')+ theme(aspect.ratio = .8)
chao
Rich = ggplot(af, aes(Age, Richness, colour = SampleType)) + geom_point()+xlab('DoL')+
  scale_color_manual(name='Body Site',values=c('Gut'=brown, 'Oral cavity'=pink, 'Skin'=orange))+
  scale_fill_manual(name='Body Site',values=c('Gut'=brown, 'Oral cavity'=pink, 'Skin'=orange))+
  geom_smooth(method='lm') + facet_wrap(~SampleType)+  stat_cor(method = 'pearson', label.x = 0.1,label.y=18) +
  theme(strip.text.x = element_text(size = 14), axis.title= element_text(size =18),axis.text = element_text(size=12),
        legend.position = 'none') + theme(aspect.ratio = .8)
Rich
#final <- grid_arrange_shared_legend(list(shan,chao,Rich))
final = grid.arrange(shan, chao, Rich, nrow=3)
ggsave(paste0(exportDir,'AlphaDiversitiesOnePlot.svg'), final, height = 8, width=10, dpi=200)
ggsave(paste0(exportDir,'AlphaDiversitiesOnePlot.png'), final, height = 8, width=10, dpi=200)

