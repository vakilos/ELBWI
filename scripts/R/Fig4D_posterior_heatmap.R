library(pheatmap)
library(ComplexHeatmap)

## work in the tmp directory
## save figures to exportDir
exportDir = ''


df = read.csv('PosteriorHeatmap.csv', header = TRUE, row.names = "ASV", check.names = FALSE )
pf = read.csv('BodysitePrevalenceAsPercentage.csv',header=TRUE, row.names="ASV", check.names = FALSE)

H = Heatmap(as.matrix(df), name="Conditional probability", show_row_names = FALSE, column_names_gp = gpar(fontsize = 10), 
        col = circlize::colorRamp2(c(0, 1), c("gold", "darkblue")), width=2,  top_annotation = HeatmapAnnotation(hist = anno_histogram(df),height=unit(2,'cm')))+ 
  HeatmapAnnotation(points = anno_points(df, which='row'), which='row') +
Heatmap(as.matrix(pf), name='Body site prevalence (%)', 
        col = circlize::colorRamp2(c(0, 1), c("lightblue", "purple")), 
        width = .5, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10))
H

pdf(paste0(exportDir,"Fig4D_posterior_heatmap.pdf"), width=11, height=8)
H
dev.off()
