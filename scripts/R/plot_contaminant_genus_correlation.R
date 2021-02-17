library(data.table)
library(ggplot2)
library(forcats)
## work in the project directory


df = fread('/media/christos/ssd/work/Infants/publication/cor_calda.csv')
ggplot(df, aes(x=fct_reorder(genus, -r),y=r, fill=contaminant)) + geom_col(width = 0.6, color='black') +
  theme(axis.text.x = element_text(angle = 90, vjust=0.0)) + scale_fill_manual(values=c('yes'='#e75353', 'no'='grey')) +
  labs(x='Genus', y='Pearson correlation coeffient (r)')
ggsave("Other_reagent_contaminants_correlate_with_Caldalkalibacillus.svg")
