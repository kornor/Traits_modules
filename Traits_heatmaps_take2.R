
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
cormat2 <- read.table("Correlation_matrix_TCGA_traits_trimmed.txt", sep = "\t",
                      header = TRUE, row.names = 1)
### for the METABRIC traits
cormat3 <- read.table("Correlation_matrix_METABRIC_traits_trimmed.txt", sep = "\t",
                      header = TRUE, row.names = 1)





############## trying out different heatmaps for the traits


#### information about par
### mai c(bottom, left, top, right) - the margin in inches
##  mar is as mai but giving the value in lines
### oma c(bottom, left, top, right) of outer margins 
par(oma= c(3,1,1,1))

### Info about heatmap.2
### must be done on a matrix not a dataframe
### dendro can be set to none or some
### key gives you a colour key
### srtCol/ srtRow gives an angle for the text
##sepCol, etc gives a cell outline
### RowSideColors - can set colours for marking the cols or rows



### colorpanel in gplots
###redgreen(n)
###redblue(n)
###bluered(n)
###colorpanel(n, low, mid, high)
## or can set the colour panel before
##eg palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)


heatmap.2(as.matrix(cormat2),dendrogram="none", 
          notecol="black",col=bluered(200),scale="none",
          key=TRUE, keysize=1.5,key.title = NA, key.xlab = "Correlation",
          density.info="none", trace="none",
          cexRow=1.2,cexCol=1.2, srtCol = 75,
          main = "Traits expression \n by module in TCGA dataset",
          rowsep = 1:nrow(cormat2), colsep = 1:ncol(cormat2), sepcolor = "grey")


heatmap.2(as.matrix(cormat3),dendrogram="none", 
          notecol="black",col=bluered(200),scale="none",
          key=TRUE, key.title = NA, key.xlab = "Correlation", keysize=1.5,
          density.info="none", trace="none",
          cexRow=1.2,cexCol=1.2, srtCol = 75,
          main = "Traits expression \n by module in METABRIC dataset",
          rowsep = 1:nrow(cormat3), colsep = 1:ncol(cormat3), sepcolor = "grey")

#################### the meth heatmap
### pheatmap function uses slightly different wording
### annotate_rows instead of 'rowsidecols' for eg



pdf("Heat_test_meth3.pdf",height=8,width=10)
fontsize = 10

pheatmap(t(cormat1), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = T, 
         treeheight_col = 0,
         main = "Methylation of Stirzaker genes by module")
dev.off()


heatmap.2(as.matrix(t(cormat1)),dendrogram="row", 
          notecol="black",col=colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100),
          scale="none",key=TRUE, keysize=1.5,
          density.info="none", trace="none",
          labRow = NA, cexRow=1.2,srtCol = 75,
          main = "Methylation of Stirzaker genes by module")

########### cormat 4 might need a different approach??

#### for the lncs vs other
cormat4 <- read.table("LncMod_cor_MetaMod.txt", sep = "\t",
                      header = TRUE, row.names = 1)
#### for the trimmed out lnc mods
cormat5 <- read.table("LncMod_cor_MetaMod_trimmed.txt", sep = "\t",
                      header = TRUE, row.names = 1)

par(oma = c(2,2,2,3))
rowstrip <- c("Turquoise", "Yellow","Green","Magenta")
colstrip <- c("Blue","Brown", "Green", "Tan", "Yellow")

heatmap.2(as.matrix(cormat5),dendrogram="none", 
          Colv = TRUE, Rowv = FALSE,
          notecol="black",col=colorpanel(200, "Green", "white", "Red"),scale="none",
          key=TRUE, key.title = NA, key.xlab = "Correlation", keysize=1.5,
          density.info="none", trace="none",
          cexRow=1 ,cexCol=1, srtCol = 75,
          main = "Correlation between \n LNC modules and EXP modules",
          rowsep = 1:nrow(cormat5), colsep = 1:ncol(cormat5), sepcolor = "grey")
