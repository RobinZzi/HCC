library(pheatmap)
library(dplyr)
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs/methy_mtx/100kbin")

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/methy_mtx")
pmd_anno <- fread("ga_global_PMD.bedGraph")
pmd_row_anno <- select(pmd_anno,V4,V9)
pmd_row_anno <- as.data.frame(pmd_row_anno[,-1])
row.names(pmd_row_anno) <- pmd_anno$V4
colnames(pmd_row_anno) <- "region_type"



hep_100k_c1 <- fread("HepG2-shNull-1_100kbin.txt")
hep_100k_c2 <- fread("HepG2-shNull-2_100kbin.txt")
hep_100k_c3 <- fread("HepG2-shNull-3_100kbin.txt")
hep_100k_e1 <- fread("HepG2-shGADD45A-1_100kbin.txt")
hep_100k_e2 <- fread("HepG2-shGADD45A-2_100kbin.txt")
hep_100k_e3 <- fread("HepG2-shGADD45A-3_100kbin.txt")

hep_100k <- cbind(hep_100k_c1,hep_100k_c2[,6],hep_100k_c3[,6],hep_100k_e1[,6],hep_100k_e2[,6],hep_100k_e3[,6])
hep_100k_sum <- as.data.frame(hep_100k[,6:11])
row.names(hep_100k_sum) <- hep_100k$id
hep_100k_sum <- na.omit(hep_100k_sum)
pheatmap(hep_100k_sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 10,
         angle_col = 45,
         border=FALSE,annotation_row = pmd_row_anno) 

hep_100k_baseMean = as.data.frame(colMeans(hep_100k_sum) ) 



snu_100k_c1 <- fread("SNU449-shNull-1_100kbin.txt")
snu_100k_c2 <- fread("SNU449-shNull-2_100kbin.txt")
snu_100k_c3 <- fread("SNU449-shNull-3_100kbin.txt")
snu_100k_e1 <- fread("SNU449-shGADD45A-1_100kbin.txt")
snu_100k_e2 <- fread("SNU449-shGADD45A-2_100kbin.txt")
snu_100k_e3 <- fread("SNU449-shGADD45A-3_100kbin.txt")

snu_100k <- cbind(snu_100k_c1,snu_100k_c2[,6],snu_100k_c3[,6],snu_100k_e1[,6],snu_100k_e2[,6],snu_100k_e3[,6])
snu_100k_sum <- as.data.frame(snu_100k[,6:11])
row.names(snu_100k_sum) <- snu_100k$id
snu_100k_sum <- na.omit(snu_100k_sum)

pheatmap(snu_100k_sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 10,
         angle_col = 45,
         border=FALSE,annotation_row = pmd_row_anno) 
  snu_100k_baseMean = as.data.frame(colMeans(snu_100k_sum) ) 




setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs/methy_mtx/10kbin")


hep_10k_c1 <- fread("HepG2-shNull-1_10kbin.txt")
hep_10k_c2 <- fread("HepG2-shNull-2_10kbin.txt")
hep_10k_c3 <- fread("HepG2-shNull-3_10kbin.txt")
hep_10k_e1 <- fread("HepG2-shGADD45A-1_10kbin.txt")
hep_10k_e2 <- fread("HepG2-shGADD45A-2_10kbin.txt")
hep_10k_e3 <- fread("HepG2-shGADD45A-3_10kbin.txt")

hep_10k <- cbind(hep_10k_c1,hep_10k_c2[,6],hep_10k_c3[,6],hep_10k_e1[,6],hep_10k_e2[,6],hep_10k_e3[,6])
hep_10k_sum <- as.data.frame(hep_10k[,6:11])
row.names(hep_10k_sum) <- hep_10k$id
hep_10k_sum <- na.omit(hep_10k_sum)
pheatmap(hep_10k_sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE) 
hep_10k_baseMean = as.data.frame(colMeans(hep_10k_sum) ) 


snu_10k_c1 <- fread("SNU449-shNull-1_10kbin.txt")
snu_10k_c2 <- fread("SNU449-shNull-2_10kbin.txt")
snu_10k_c3 <- fread("SNU449-shNull-3_10kbin.txt")
snu_10k_e1 <- fread("SNU449-shGADD45A-1_10kbin.txt")
snu_10k_e2 <- fread("SNU449-shGADD45A-2_10kbin.txt")
snu_10k_e3 <- fread("SNU449-shGADD45A-3_10kbin.txt")

snu_10k <- cbind(snu_10k_c1,snu_10k_c2[,6],snu_10k_c3[,6],snu_10k_e1[,6],snu_10k_e2[,6],snu_10k_e3[,6])
snu_10k_sum <- as.data.frame(snu_10k[,6:11])
row.names(snu_10k_sum) <- snu_10k$id
snu_10k_sum <- na.omit(snu_10k_sum)
pheatmap(snu_10k_sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE) 
snu_10k_baseMean = as.data.frame(colMeans(snu_10k_sum) ) 




sum <- cbind(snu_100k_c1,snu_100k_c2[,6],snu_100k_c3[,6],snu_100k_e1[,6],snu_100k_e2[,6],snu_100k_e3[,6],
             hep_100k_c1[,6],hep_100k_c2[,6],hep_100k_c3[,6],hep_100k_e1[,6],hep_100k_e2[,6],hep_100k_e3[,6])
sum <- sum[,6:17]
row.names(sum) <- hep_100k$id
sum <- na.omit(sum)
pheatmap(sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE,annotation_row = pmd_row_anno) 




pmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "PMD"))
HEPPMDMean = as.data.frame(colMeans(hep_100k_sum[pmd_regions,]) ) 
colnames(HEPPMDMean) <- "HepG2-PMD"
SNUPMDMean = as.data.frame(colMeans(snu_100k_sum[pmd_regions,]) ) 
colnames(SNUPMDMean) <- "SNU449-PMD"

colnames(G6baseMean) <- "G6-base"
colnames(GAbaseMean) <- "GA-base"








pmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "PMD"))
hmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "HMD"))
snu_100k_sum <-na.omit(snu_100k_sum)
HepPMDMean = as.data.frame(colMeans(hep_100k_sum[pmd_regions,]) ) 
colnames(HepPMDMean) <- "Hep-PMD"
SnuPMDMean = as.data.frame(colMeans(snu_100k_sum[hmd_regions,],na.rm = T) ) 
colnames(SnuPMDMean) <- "Snu-PMD"
HepHMDMean = as.data.frame(colMeans(hep_100k_sum[hmd_regions,]) ) 
colnames(HepHMDMean) <- "Hep-HMD"
SnuHMDMean = as.data.frame(colMeans(snu_100k_sum[hmd_regions,],na.rm = T) ) 
colnames(SnuHMDMean) <- "Snu-HMD"

colnames(snu_100k_baseMean) <- "Snu-global"
colnames(hep_100k_baseMean) <- "Hep-global"

Hepmeansum <- cbind(HepPMDMean,HepHMDMean,hep_100k_baseMean)
Snumeansum <- cbind(SnuPMDMean,SnuHMDMean,snu_100k_baseMean)
