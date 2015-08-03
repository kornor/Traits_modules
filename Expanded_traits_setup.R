###

######  THIS IS A SCRIPT TO SET UP TRAITS AGAINST THE TCGA AND META WGCNA MODULES

#set the working directory and open the libraries needed

setwd("~/Bioinformatics Work/Meth & RNA/Traits_Modules_Work")

library(plyr)
library(survival)
library(ggplot2)
library(gplots)
library(WGCNA)
library(pheatmap)
library(RColorBrewer)

#### load in the files needed

load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")

###create a table of the module colour by gene 
colorsA1 = names(table(modules1))



##### load in the traits list of interest
## Pam50

pam50 <- read.table("Pam50_geneList.txt", sep = "\t", header = TRUE)

##  Subset the TCGA set with these genes

list <- intersect(rownames(datExpr1), pam50[,1])
pam_t <- datExpr1[list,]
pam_t <- t(pam_t)
### and the METABRIC

list <- intersect(rownames(datExpr2), pam50[,1])
pam_m <- datExpr2[list,]
pam_m <- t(pam_m)

### Burstein set
  bur80 <- read.table("Burstein_TNBC_signature.txt", sep = "\t", header = TRUE)

##  Subset the TCGA set with these genes

list <- intersect(rownames(datExpr1), bur80[,1])
bur_t <- datExpr1[list,]
  bur_t <- t(bur_t)
### and the METABRIC

list <- intersect(rownames(datExpr2), bur80[,1])
bur_m <- datExpr2[list,]
  bur_m <- t(bur_m)
#################################

##### ok what now?  
#####  Do correlations between module expression and the sets??

### try with the TCGA pam50
  nGenes = ncol(datExpr1);
  nSamples = nrow(datExpr1);
moduleTraitCor1 = cor(ME_1A, pam_t, use = "p");
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamples);

  
## recall that setting treeheight_col or _row to '0' will remove 
  ##dendro for that side
  pheatmap(t(moduleTraitCor1), 
           color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
           kmeans_k = NA, breaks = NA, border_color = "grey60",
           cellwidth = NA, cellheight = NA, scale = "none", 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           fontsize_row = 10,
           annotation_legend = TRUE, drop_levels = TRUE, 
           show_rownames = F,show_colnames = T, 
           main = "Expression of Pam50 genes by module (TCGA)")

  ######################  
  ### try with the METABRIC pam50
  nGenes = ncol(datExpr2);
  nSamples = nrow(datExpr2);
  moduleTraitCor2 = cor(ME_2A, pam_m, use = "p");
  moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamples);
  
  
  ## recall that setting treeheight_col or _row to '0' will remove 
  ##dendro for that side
  pheatmap(t(moduleTraitCor2), 
           color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
           kmeans_k = NA, breaks = NA, border_color = "grey60",
           cellwidth = NA, cellheight = NA, scale = "none", 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           fontsize_row = 10,
           annotation_legend = TRUE, drop_levels = TRUE, 
           show_rownames = F,show_colnames = T, 
           main = "Expression of Pam50 genes by module (METABRIC)")
  
  
  
############# Bur set
  
  
  ### try with the TCGA pam50
  nGenes = ncol(datExpr1);
  nSamples = nrow(datExpr1);
  moduleTraitCor3 = cor(ME_1A, bur_t, use = "p");
  moduleTraitPvalue3 = corPvalueStudent(moduleTraitCor3, nSamples);
  
  
  pdf("Heat_bur80.pdf",height=8,width=10)
  fontsize = 10
  
  ## recall that setting treeheight_col or _row to '0' will remove 
  ##dendro for that side
  pheatmap(t(moduleTraitCor3), 
           color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
           kmeans_k = NA, breaks = NA, border_color = "grey60",
           cellwidth = NA, cellheight = NA, scale = "none", 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           fontsize_row = 5,
           annotation_legend = TRUE, drop_levels = TRUE, 
           show_rownames = T,show_colnames = T, 
           main = "Expression of Bur80 genes by module (TCGA)")
  
  ######################  
  ### try with the METABRIC bur
  
  
  nGenes = ncol(datExpr2);
  nSamples = nrow(datExpr2);
  moduleTraitCor4 = cor(ME_2A, bur_m, use = "p");
  moduleTraitPvalue4 = corPvalueStudent(moduleTraitCor4, nSamples);
  
  

  
  ## recall that setting treeheight_col or _row to '0' will remove 
  ##dendro for that side
  pheatmap(t(moduleTraitCor4), 
           color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
           kmeans_k = NA, breaks = NA, border_color = "grey60",
           cellwidth = NA, cellheight = NA, scale = "none", 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           fontsize_row = 5,
           annotation_legend = TRUE, drop_levels = TRUE, 
           show_rownames = T,show_colnames = T, 
           main = "Expression of Bur80 genes by module (METABRIC)")
  dev.off()
  
  
  
#### load the individual Bur type signatures
  
  bur4 <- read.table("Bur80_4_subtypes.txt", sep = "\t", header = TRUE)
### LAR
  list <- intersect(rownames(datExpr2), bur4[,1])
  lar_m <- datExpr2[list,]
### MES  
  list <- intersect(rownames(datExpr2), bur4[,2])
  mes_m <- datExpr2[list,]
## BLIA
  list <- intersect(rownames(datExpr2), bur4[,3])
  blia_m <- datExpr2[list,]
## BLIS
  list <- intersect(rownames(datExpr2), bur4[,4])
  blis_m <- datExpr2[list,]
  
####get the mean centred ave of all of these sample
  # create function to center : 'colMeans()'
  center_colmeans <- function(x) {
    xcenter = colMeans(x)
    x - rep(xcenter, rep.int(nrow(x), ncol(x)))
  }
  
  centre_t <- center_colmeans(t(lar_m))
  
  meta_ave <- as.data.frame(rowMeans(centre_t))
  
  rownames(meta_ave) <- colnames(lar_m) 
  colnames(meta_ave)[1] <- "LAR signature"
  
  
  ### for MES
  
  centre_t <- center_colmeans(t(mes_m))
  meta_ave$MES_signature <- rowMeans(centre_t)
  
  ## for BLIA
  
  centre_t <- center_colmeans(t(blia_m))
  meta_ave$BLIA_signature <- rowMeans(centre_t)
  
  ## for BLIS
  
  centre_t <- center_colmeans(t(blis_m))
  meta_ave$BLIS_signature <- rowMeans(centre_t)
  
  
##### META_AVE is now the mean-centred average expression
  ## of each signautre for each sample
  
  ##Try correlations again
  
  ### try with the TCGA pam50
  nGenes = ncol(ME_2A);
  nSamples = nrow(ME_2A);
  moduleTraitCor5 = cor(ME_2A, meta_ave, use = "p");
  moduleTraitPvalue5 = corPvalueStudent(moduleTraitCor5, nSamples);
  
  
  pdf("Heat_bur_subtype_meta.pdf",height=8,width=10)
  fontsize = 10
  
  ## recall that setting treeheight_col or _row to '0' will remove 
  ##dendro for that side
  pheatmap(t(moduleTraitCor5), 
           color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
           kmeans_k = NA, breaks = NA, border_color = "grey60",
           cellwidth = NA, cellheight = NA, scale = "none", 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           fontsize_row = 10,
           annotation_legend = TRUE, drop_levels = TRUE, 
           show_rownames = T,show_colnames = T, 
           main = "Expression of Bur80 genes by module (TCGA)")
  
 dev.off() 
 
 write.table(meta_ave, "Burstein_METABRIC.txt", sep = "\t")
 