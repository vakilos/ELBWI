library(ggplot2)
library(data.table)
library(ggsignif)
library(ggpval)
library(ggpubr)

## work in the tmp directory
## save figures to exportDir

df = fread('DistanceBodysites_toR.tsv')
exportDir = ''

### wilcoxon rank sum test - same with Mann Whitney U test (two-sided (specify it in python))
compare_means(Distance ~ Comparison, group.by ="Patient", data = df) 

plt = ggplot(df, aes(x=Patient, y=Distance)) + geom_boxplot(aes(fill=Comparison),position='dodge', outlier.size = .2) + 
theme(legend.position = 'bottom' ,legend.text=element_text(size=14),aspect.ratio = 0.5 ,
      axis.title.y = element_text(size = 18),axis.title.x = element_text(size = 18), axis.text = element_text(size = 12) ) + 
ylab('Bray-Curtis distance') +xlab("Infant") + scale_fill_discrete("",labels= c("between body sites", 'within body site')) +
geom_text(size=6,data = data.frame(Patient = c('P01',"P02",'P03','P15'), Distance=c(1.05,1.05,1.05,1.05)), label= c("*","***","*","***") )
plt



ggsave(paste0(exportDir,'BodysiteDistance.svg'),dpi=200, height=9, width = 10)
ggsave(paste0(exportDir,'BodysiteDistance.png'),dpi=200, height=9, width = 10)




