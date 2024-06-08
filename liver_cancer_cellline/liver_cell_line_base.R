
setwd("~/projects/hcc/scripts/liver_cancer_cellline")
rm(list=ls())
save.image("diff.Rdata")
mat = read.delim("combined-chrM.counts",header=T,skip=1)
counts = mat[,-c(1:6)]
rownames(counts) = mat$Geneid
list <- c("MET","GADD45A","SNHG6","ONECUT2")

list_mat <- counts[list,]



rpkm_mat = read.delim("combined-chrM.rpkm",header=F,skip=1)
colnames(rpkm_mat) <- colnames(mat)
rpkm = rpkm_mat[,-c(1:6)]
rownames(rpkm) = rpkm_mat$Geneid
list <- c("MET","GADD45A","SNHG6","ONECUT2")

list_rpkm <- rpkm[list,]
