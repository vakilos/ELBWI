library(ggplot2)
library(data.table)
library(gtable)
library(grid)
library(svglite)

## work in the tmp directory
## save figures to exportDir

## First, create a tsv file withi envfit.R results
df = fread('envFit_to_R.tsv')
exportDir =''
df$r2 = round(df$r2,3)
MaxLimit = max(df$r2)+ 0.1*max(df$r2)
ef = subset(df, df$Table == 'exp')
pf = subset(df, df$Table == 'pre')

ef$Variable = factor(ef$Variable, levels = unique(ef$Variable[order(ef$grouping)]), ordered=TRUE)
cols = c('#21cc8e' = '#21cc8e', '#3581bd' = '#3581bd','#e01f93'='#e01f93')
p1 = ggplot(ef, aes(y=r2, x=Variable, fill=Color, label=Asterisk), color='black') + geom_bar(stat='identity', width=.5) + coord_flip() + 
  facet_wrap(~SampleType) + scale_fill_manual(values = cols) + 
  theme(strip.text.x = element_text(size = 14), axis.text=element_text(size=14), axis.title=element_text(size=16),aspect.ratio = 1,legend.position = 'none', axis.text.x = element_blank(),axis.ticks.x = element_blank()) + 
  labs(x='',y='') + geom_text(size=6,hjust=-1, vjust=0.7) +
  scale_y_continuous(limits= c(0,MaxLimit), 
                     breaks=c(round(MaxLimit/4,2),round(2*MaxLimit/4,2),round(3*MaxLimit/4,2),round(MaxLimit,2) ))
p1

ggsave(paste0(exportDir,'envfit_pt1.svg'), height=8,width=12)

p2 = ggplot(pf, aes(y=r2, x=Variable, fill='#d16666',label=Asterisk),colour='black') + geom_bar(stat = 'identity', width=.5)+ coord_flip()+
  facet_wrap(~SampleType)+ 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16), aspect.ratio = 0.6, legend.position = 'none', strip.background = element_blank(), strip.text =element_blank())+
  labs(x='',y=expression(paste('r')^2)) + geom_text(size=6,hjust=-.6, vjust=0.7) +
  scale_y_continuous(limits= c(0,MaxLimit), 
                     breaks=c(round(MaxLimit/4,2),round(2*MaxLimit/4,2),round(3*MaxLimit/4,2),round(MaxLimit,2) ))
p2

ggsave(paste0(exportDir,'envfit_pt2.svg'), height=8,width=12)

###### plot everything in one plot
g1 = ggplotGrob(p1)
g2 = ggplotGrob(p2)
g = rbind(g1,g2, size='first')
g$widths = unit.pmax(g1$widths, g2$widths) + unit.pmax(g1$widths, g2$widths)*0.05
grid.newpage()
svglite(paste0(exportDir,'envfit_final_plot.svg'), height =6, width = 10)
grid.draw(g)
dev.off()
  
