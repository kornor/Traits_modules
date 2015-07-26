
##### This is a script for all the heatmaps wahahahah!

### differnt working directory for ease

setwd("~/Bioinformatics Work/Meth & RNA/Traits_Modules_Work")

### Load libraries

library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(gplots)


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

cormat1 <- read.table("Correlation_matrix_TCGA_meth.txt", sep = "\t", header = TRUE, row.names = 1)  
cormat1 <- t(cormat1)  
  
  fontsize = 10
  
pdf("Heat_test_meth.pdf",height=8,width=10)
pheatmap(cormat1, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = TRUE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, show_rownames = F,
         show_colnames = T, main = "Methylation of genes by module")
dev.off()



pdf("Heat_test_meth1.pdf",height=8,width=10)
pheatmap(t(cormat1), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = TRUE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T,
         show_colnames = F, main = "Methylation of genes by module")
dev.off()


pdf("Heat_test2.pdf",height=8,width=8)
pheatmap(cormat2, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = T,
         cluster_cols = FALSE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T,
         show_colnames = T, main = "Correlation of expression of \n trait genes by module (TCGA)")
dev.off()

pdf("Heat_test3.pdf",height=8,width=8)
pheatmap(cormat3, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = T,
         cluster_cols = FALSE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T,
         show_colnames = T, main = "Correlation of expression of \n trait genes by module (METABRIC)")
dev.off()

