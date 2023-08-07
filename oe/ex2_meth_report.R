
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex2/methy_report")

rm(list=ls())
library(data.table)
library(R.utils)
library(tidyr)
library(pheatmap)
library(edgeR)
library(stringr)
library(FactoMineR)
library(dplyr)


load("methy_merge.Rdata")
save.image("methy_merge.Rdata")

save.image("methy_merge.Rdata")
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
          "chr16","chr17","chr18","chr19","chr20","chr21","chr22")
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex2/methy_report/report")
fs <- list.files()
for (i in 1:length(fs)) {
  setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex2/methy_report/report")
  data <- fread(paste(fs[i]))
  data <- subset(data,data$V1 %in% chrs)
  colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
  data[,3] <- round(data[,2]/100000)
  data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
  data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
  sample_id <-str_split_fixed(fs[i],"_depulicated.CpG_report.txt",2)[1]
  data[,paste(sample_id)] <- data[,2]/(data[,2]+data[,3])
  setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex2/methy_report")
  write.table(data,paste(sample_id,".txt",sep = ""))
}


GA_C1 <- fread("HepG2-GADD45A-C1.txt")
GA_C2 <- fread("HepG2-GADD45A-C2.txt")
GA_C3 <- fread("HepG2-GADD45A-C3.txt")
GA_E1 <- fread("HepG2-GADD45A-OE1.txt")
GA_E2 <- fread("HepG2-GADD45A-OE2.txt")
GA_E3 <- fread("HepG2-GADD45A-OE3.txt")

GA <- cbind(GA_C1,GA_C2[,5],GA_C3[,5],GA_E1[,5],GA_E2[,5],GA_E3[,5])
GA_mtx <- as.data.frame(GA[,5:10])
row.names(GA_mtx) <- GA$id
GA_mtx <- na.omit(GA_mtx)
pheatmap(GA_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE) 
GAbaseMean = as.data.frame(colMeans(GA_mtx) ) 


G6_C1 <- fread("HepG2-SNHG6-C1.txt")
G6_C2 <- fread("HepG2-SNHG6-C2.txt")
G6_C3 <- fread("HepG2-SNHG6-C3.txt")
G6_E1 <- fread("HepG2-SNHG6-OE1.txt")
G6_E2 <- fread("HepG2-SNHG6-OE2.txt")
G6_E3 <- fread("HepG2-SNHG6-OE3.txt")

G6 <- cbind(G6_C1,G6_C2[,5],G6_C3[,5],G6_E1[,5],G6_E2[,5],G6_E3[,5])
G6_mtx <- as.data.frame(G6[,5:10])
row.names(G6_mtx) <- G6$id
G6_mtx <- na.omit(G6_mtx)
pheatmap(G6_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,annotation_row = pmd_row_anno,
         angle_col = 45,treeheight_row = 0,
         border=FALSE) 
G6baseMean = as.data.frame(colMeans(G6_mtx) )














