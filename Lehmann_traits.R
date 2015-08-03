#### script for Lehmann subtype comparison




#### load the individual Bur type signatures

leh6 <- read.table("Lehmann_subtype_genes.txt", sep = "\t", header = TRUE)
### BL1
list <- intersect(rownames(datExpr2), leh6[,1])
bl1_m <- datExpr2[list,]
### BL2  
list <- intersect(rownames(datExpr2), leh6[,2])
bl2_m <- datExpr2[list,]
## IM
list <- intersect(rownames(datExpr2), leh6[,3])
im_m <- datExpr2[list,]
## M
list <- intersect(rownames(datExpr2), leh6[,4])
m_m <- datExpr2[list,]
### MSL
list <- intersect(rownames(datExpr2), leh6[,5])
msl_m <- datExpr2[list,]
### LAR
list <- intersect(rownames(datExpr2), leh6[,6])
lar_m <- datExpr2[list,]

##################################################
####get the mean centred ave of all of these sample
# create function to center : 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

centre_t <- center_colmeans(t(lar_m))

meta_ave <- as.data.frame(rowMeans(centre_t))

rownames(meta_ave) <- colnames(lar_m) 
colnames(meta_ave)[1] <- "LAR_signature"


### for M

centre_t <- center_colmeans(t(m_m))
meta_ave$M_signature <- rowMeans(centre_t)

## for MSL

centre_t <- center_colmeans(t(msl_m))
meta_ave$MSL_signature <- rowMeans(centre_t)

## for IM

centre_t <- center_colmeans(t(im_m))
meta_ave$IM_signature <- rowMeans(centre_t)

## for BL1

centre_t <- center_colmeans(t(bl1_m))
meta_ave$BL1_signature <- rowMeans(centre_t)

## for BL2

centre_t <- center_colmeans(t(bl2_m))
meta_ave$BL2_signature <- rowMeans(centre_t)


##### META_AVE is now the mean-centred average expression
## of each signautre for each sample
## write out

write.table(meta_ave, "Lehmann_METABRIC.txt", sep = "\t")
##Try correlations again

### try with the TCGA pam50
nGenes = ncol(ME_2A);
nSamples = nrow(ME_2A);
moduleTraitCor5 = cor(ME_2A, meta_ave, use = "p");
moduleTraitPvalue5 = corPvalueStudent(moduleTraitCor5, nSamples);


pdf("Heat_leh_subtype_meta.pdf",height=8,width=10)
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
         main = "Expression of Lehmann subtype by module (METABRIC)")

dev.off() 
