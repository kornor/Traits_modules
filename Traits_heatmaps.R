
##### This is a script for all the heatmaps wahahahah!

### differnt working directory for ease

setwd("~/Bioinformatics Work/Meth & RNA/Traits_Modules_Work")

### Load libraries

library(RColorBrewer)
library(pheatmap)
library(ggplot2)



### need a cormat for the modules vs Stir meth

  cormat1 <- read.table("Correlation_matrix_TCGA_meth.txt", sep = "\t",
                      header = TRUE, row.names = 1)
### for the TCGA and traits
  cormat2 <- read.table("Correlation_matrix_TCGA_traits.txt", sep = "\t",
                      header = TRUE, row.names = 1)
### for the METABRIC traits
  cormat3 <- read.table("Correlation_matrix_METABRIC_traits.txt", sep = "\t",
                      header = TRUE, row.names = 1)



## want a heat map of the methylation v modules that shows the ordering etc
### not sure about that for the other traits?

  fontsize = 10
  
pheatmap(cormat1, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T,
         show_colnames = F, main = "Methylation of genes by module")

