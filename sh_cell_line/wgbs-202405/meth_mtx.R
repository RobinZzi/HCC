library(pheatmap)
library(dplyr)
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/methy_mtx/100kbin")

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/methy_mtx")
pmd_anno <- fread("ga_global_PMD.bedGraph")
pmd_row_anno <- select(pmd_anno,V4,V9)
pmd_row_anno <- as.data.frame(pmd_row_anno[,-1])
row.names(pmd_row_anno) <- pmd_anno$V4
colnames(pmd_row_anno) <- "region_type"



huh_100k_e1 <- fread("Huh7_SH2-1_100kbin.txt")
huh_100k_e2 <- fread("Huh7_SH2-2_100kbin.txt")
huh_100k_c1 <- fread("Huh7_V1_100kbin.txt")
huh_100k_c2 <- fread("Huh7_V2_100kbin.txt")


huh_100k <- cbind(huh_100k_e1,huh_100k_e2[,6],huh_100k_c1[,6],huh_100k_c2[,6])
huh_100k_sum <- as.data.frame(huh_100k[,6:9])
row.names(huh_100k_sum) <- huh_100k$id
huh_100k_sum <- na.omit(huh_100k_sum)
pheatmap(huh_100k_sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 10,
         angle_col = 45,
         border=FALSE,annotation_row = pmd_row_anno) 

huh_100k_baseMean = as.data.frame(colMeans(huh_100k_sum) ) 









setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/methy_mtx/10kbin")


huh_10k_e1 <- fread("Huh7_SH2-1_10kbin.txt")
huh_10k_e2 <- fread("Huh7_SH2-2_10kbin.txt")
huh_10k_c1 <- fread("Huh7_V1_10kbin.txt")
huh_10k_c2 <- fread("Huh7_V2_10kbin.txt")

huh_10k <- cbind(huh_10k_e1,huh_10k_e2[,6],huh_10k_c1[,6],huh_10k_c2[,6])
huh_10k_sum <- as.data.frame(huh_10k[,6:9])
row.names(huh_10k_sum) <- huh_10k$id
huh_10k_sum <- na.omit(huh_10k_sum)
pheatmap(huh_10k_sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 10,
         angle_col = 45,
         border=FALSE,annotation_row = pmd_row_anno) 

huh_10k_baseMean = as.data.frame(colMeans(huh_100k_sum) ) 





pmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "PMD"))
huh_100k_pmd <- huh_100k_sum[pmd_regions,]
huh_100k_pmd <- na.omit(huh_100k_pmd)
huh_100k_PMDMean = as.data.frame(colMeans(huh_100k_pmd) ) 
colnames(Huh_100k_PMDMean) <- "Huh7-PMD"


colnames(G6baseMean) <- "G6-base"
colnames(GAbaseMean) <- "GA-base"