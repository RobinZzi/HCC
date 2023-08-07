
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/methy_merge")

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
save.image("methy_merge_ex1.Rdata")

save.image("methy_merge.Rdata")
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16","chr17","chr18","chr19","chr20","chr21","chr22")

fs <- list.files()
for (i in 1:length(fs)) {
setwd("~/projects/hcc/analysis/oe/meth_ex/methy_report")
data <- fread(paste(fs[i]))
data <- subset(data,data$V1 %in% chrs)
colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
data[,3] <- round(data[,2]/100000)
data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
sample_id <-str_split_fixed(fs[i],"_remove_dup.CpG_report.txt",2)[1]
data[,paste(sample_id)] <- data[,2]/(data[,2]+data[,3])
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/methy_merge")
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





plot(1:10)









G6_t <- t(G6_mtx)
G6_group <- data.frame(Sample = rownames(G6_t), Group = rep(c("Control", "OE"), each = 3))
G6_pca <- PCA(G6_t, ncp = 2, scale.unit = TRUE, graph = FALSE)
G6_pca_sample <- data.frame(G6_pca$ind$coord[ ,1:2])
G6_pca_sample$Sample=row.names(G6_pca_sample)
G6_pca_eig1 <- round(G6_pca$eig[1,2], 2)
G6_pca_eig2 <- round(G6_pca$eig[2,2],2 )


G6_pca_sample <- merge(G6_pca_sample,G6_group,by="Sample")
head(G6_pca_sample)
G6_p <- ggplot(data = G6_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 2) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', G6_pca_eig1, '%'), y = paste('PCA2:', G6_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="SNHG6_pca")


GA_t <- t(GA_mtx)
GA_group <- data.frame(Sample = rownames(GA_t), Group = rep(c("Control", "OE"), each = 3))
GA_pca <- PCA(GA_t, ncp = 2, scale.unit = TRUE, graph = FALSE)
GA_pca_sample <- data.frame(GA_pca$ind$coord[ ,1:2])
GA_pca_sample$Sample=row.names(GA_pca_sample)
GA_pca_eig1 <- round(GA_pca$eig[1,2], 2)
GA_pca_eig2 <- round(GA_pca$eig[2,2],2 )


GA_pca_sample <- merge(GA_pca_sample,GA_group,by="Sample")
head(GA_pca_sample)
GA_p <- ggplot(data = GA_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 2) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', GA_pca_eig1, '%'), y = paste('PCA2:', GA_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="GADD45A_pca")








row.names(GA_mtx) <- GA$id



global_bed <- as.data.frame(str_split_fixed(row.names(GA_mtx),"_",2))
colnames(global_bed) <- c("chr","start")
global_bed$cite <- row.names(GA_mtx)
global_bed$start <- as.numeric(global_bed$start)
global_bed$end <- global_bed$start+1
global_bed$start <- global_bed$start*100000
global_bed$end <- global_bed$end*100000
global_bed$end <- as.numeric(global_bed$end)
global_bed$start <- as.numeric(global_bed$start)
global_bed <- select(global_bed,chr,start,end)
write.table(global_bed, "global.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


##t-test
G6_p_values <- apply(G6_mtx, 1, function(row) {
  t.test(row[1:3], row[4:6])$p.value
})
G6_fc <- apply(G6_mtx, 1, function(row) {
  mean(row[1:3])/mean(row[4:6])
})

G6_adjusted_p_values <- p.adjust(G6_p_values, method  = "fdr")
G6_adjusted_p_values_b <- p.adjust(G6_p_values, method  = "bonferroni")
G6_test <- cbind(G6_p_values,G6_adjusted_p_values,G6_adjusted_p_values_b,G6_fc)
G6_test <- as.data.frame(G6_test)
colnames(G6_test) <- c("pvalue","bf_adj","fdr","fc")
G6_test_significant <- subset(G6_test,subset =  pvalue< 0.05)
##t-test
G6_p_values <- apply(G6_mtx, 1, function(row) {
  t.test(row[1:3], row[4:6])$p.value
})
G6_fc <- apply(G6_mtx, 1, function(row) {
  mean(row[1:3])/mean(row[4:6])
})

G6_adjusted_p_values <- p.adjust(G6_p_values, method  = "fdr")
G6_adjusted_p_values_b <- p.adjust(G6_p_values, method  = "bonferroni")
G6_test <- cbind(G6_p_values,G6_adjusted_p_values,G6_adjusted_p_values_b,G6_fc)
G6_test <- as.data.frame(G6_test)
colnames(G6_test) <- c("pvalue","bf_adj","fdr","fc")
G6_test_significant <- subset(G6_test,subset =  pvalue< 0.05)
GA_p_values <- apply(GA_mtx, 1, function(row) {
  t.test(row[1:3], row[4:6])$p.value
})
GA_fc <- apply(GA_mtx, 1, function(row) {
  mean(row[1:3])/mean(row[4:6])
})

GA_adjusted_p_values <- p.adjust(GA_p_values, method  = "fdr")
GA_adjusted_p_values_b <- p.adjust(GA_p_values, method  = "bonferroni")
GA_test <- cbind(GA_p_values,GA_adjusted_p_values,GA_adjusted_p_values_b,GA_fc)
GA_test <- as.data.frame(GA_test)
colnames(GA_test) <- c("pvalue","bf_adj","fdr","fc")
GA_test_significant <- subset(GA_test,subset =  pvalue< 0.05)

GA_significant_bed <- as.data.frame(str_split_fixed(row.names(GA_test_significant),"_",2))
colnames(GA_significant_bed) <- c("chr","start")
GA_significant_bed$cite <- row.names(GA_test_significant)
GA_significant_bed$start <- as.numeric(GA_significant_bed$start)
GA_significant_bed$end <- GA_significant_bed$start+1
GA_significant_bed$start <- GA_significant_bed$start*100000
GA_significant_bed$end <- GA_significant_bed$end*100000
GA_significant_bed$end <- as.numeric(GA_significant_bed$end)
GA_significant_bed$start <- as.numeric(GA_significant_bed$start)
GA_significant_bed <- select(GA_significant_bed,chr,start,end)
write.table(GA_significant_bed, "GA_significant.bed",sep = "\t",quote=F,row.names = F,col.names = F)
G6_significant_bed <- as.data.frame(str_split_fixed(row.names(G6_test_significant),"_",2))
colnames(G6_significant_bed) <- c("chr","start")
G6_significant_bed$cite <- row.names(G6_test_significant)
G6_significant_bed$start <- as.numeric(G6_significant_bed$start)
G6_significant_bed$end <- G6_significant_bed$start+1
G6_significant_bed$start <- G6_significant_bed$start*100000
G6_significant_bed$end <- G6_significant_bed$end*100000
G6_significant_bed$end <- as.numeric(G6_significant_bed$end)
G6_significant_bed$start <- as.numeric(G6_significant_bed$start)
G6_significant_bed <- select(G6_significant_bed,chr,start,end)
write.table(G6_significant_bed, "G6_significant.bed",sep = "\t",quote=F,row.names = F,col.names = F)



pmd_anno <- fread("global_pmd.bedGraph")
pmd_row_anno <- select(pmd_anno,V4,V9)
pmd_row_anno <- as.data.frame(pmd_row_anno[,-1])
row.names(pmd_row_anno) <- pmd_anno$V4
colnames(pmd_row_anno) <- "region_type"


pheatmap(GA_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,annotation_row = pmd_row_anno,
         angle_col = 45,
         border=FALSE)



pheatmap(G6_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,annotation_row = pmd_row_anno,
         angle_col = 45,
         border=FALSE)












G6_C1_bed <- as.data.frame(str_split_fixed(G6_C1$id,"_",2))
colnames(G6_C1_bed) <- c("chr","start")
G6_C1_bed$cite <- G6_C1$id
G6_C1_bed$meth_level <- G6_C1$`HepG2-SNHG6-C1`
G6_C1_bed$start <- as.numeric(G6_C1_bed$start)
G6_C1_bed$end <- G6_C1_bed$start+1
G6_C1_bed$start <- G6_C1_bed$start*100000
G6_C1_bed$end <- G6_C1_bed$end*100000
G6_C1_bed$end <- as.numeric(G6_C1_bed$end)
G6_C1_bed$start <- as.numeric(G6_C1_bed$start)
G6_C1_bed <- select(G6_C1_bed,chr,start,end,meth_level)
write.table(G6_C1_bed, "G6_C1.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_C2_bed <- as.data.frame(str_split_fixed(G6_C2$id,"_",2))
colnames(G6_C2_bed) <- c("chr","start")
G6_C2_bed$cite <- G6_C2$id
G6_C2_bed$meth_level <- G6_C2$`HepG2-SNHG6-C2`
G6_C2_bed$start <- as.numeric(G6_C2_bed$start)
G6_C2_bed$end <- G6_C2_bed$start+1
G6_C2_bed$start <- G6_C2_bed$start*100000
G6_C2_bed$end <- G6_C2_bed$end*100000
G6_C2_bed$end <- as.numeric(G6_C2_bed$end)
G6_C2_bed$start <- as.numeric(G6_C2_bed$start)
G6_C2_bed <- select(G6_C2_bed,chr,start,end,meth_level)
write.table(G6_C2_bed, "G6_C2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_C3_bed <- as.data.frame(str_split_fixed(G6_C3$id,"_",2))
colnames(G6_C3_bed) <- c("chr","start")
G6_C3_bed$cite <- G6_C3$id
G6_C3_bed$meth_level <- G6_C3$`HepG2-SNHG6-C3`
G6_C3_bed$start <- as.numeric(G6_C3_bed$start)
G6_C3_bed$end <- G6_C3_bed$start+1
G6_C3_bed$start <- G6_C3_bed$start*100000
G6_C3_bed$end <- G6_C3_bed$end*100000
G6_C3_bed$end <- as.numeric(G6_C3_bed$end)
G6_C3_bed$start <- as.numeric(G6_C3_bed$start)
G6_C3_bed <- select(G6_C3_bed,chr,start,end,meth_level)
write.table(G6_C3_bed, "G6_C3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_E1_bed <- as.data.frame(str_split_fixed(G6_E1$id,"_",2))
colnames(G6_E1_bed) <- c("chr","start")
G6_E1_bed$cite <- G6_E1$id
G6_E1_bed$meth_level <- G6_E1$`HepG2-SNHG6-OE1`
G6_E1_bed$start <- as.numeric(G6_E1_bed$start)
G6_E1_bed$end <- G6_E1_bed$start+1
G6_E1_bed$start <- G6_E1_bed$start*100000
G6_E1_bed$end <- G6_E1_bed$end*100000
G6_E1_bed$end <- as.numeric(G6_E1_bed$end)
G6_E1_bed$start <- as.numeric(G6_E1_bed$start)
G6_E1_bed <- select(G6_E1_bed,chr,start,end,meth_level)
write.table(G6_E1_bed, "G6_E1.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_E2_bed <- as.data.frame(str_split_fixed(G6_E2$id,"_",2))
colnames(G6_E2_bed) <- c("chr","start")
G6_E2_bed$cite <- G6_E2$id
G6_E2_bed$meth_level <- G6_E2$`HepG2-SNHG6-OE2`
G6_E2_bed$start <- as.numeric(G6_E2_bed$start)
G6_E2_bed$end <- G6_E2_bed$start+1
G6_E2_bed$start <- G6_E2_bed$start*100000
G6_E2_bed$end <- G6_E2_bed$end*100000
G6_E2_bed$end <- as.numeric(G6_E2_bed$end)
G6_E2_bed$start <- as.numeric(G6_E2_bed$start)
G6_E2_bed <- select(G6_E2_bed,chr,start,end,meth_level)
write.table(G6_E2_bed, "G6_E2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_E3_bed <- as.data.frame(str_split_fixed(G6_E3$id,"_",2))
colnames(G6_E3_bed) <- c("chr","start")
G6_E3_bed$cite <- G6_E3$id
G6_E3_bed$meth_level <- G6_E3$`HepG2-SNHG6-OE3`
G6_E3_bed$start <- as.numeric(G6_E3_bed$start)
G6_E3_bed$end <- G6_E3_bed$start+1
G6_E3_bed$start <- G6_E3_bed$start*100000
G6_E3_bed$end <- G6_E3_bed$end*100000
G6_E3_bed$end <- as.numeric(G6_E3_bed$end)
G6_E3_bed$start <- as.numeric(G6_E3_bed$start)
G6_E3_bed <- select(G6_E3_bed,chr,start,end,meth_level)
write.table(G6_E3_bed, "G6_E3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



G6_C1_bed <- as.data.frame(str_split_fixed(G6_C1$id,"_",2))
colnames(G6_C1_bed) <- c("chr","start")
G6_C1_bed$cite <- G6_C1$id
G6_C1_bed$meth_level <- G6_C1$`HepG2-SNHG6-C1`
G6_C1_bed$start <- as.numeric(G6_C1_bed$start)
G6_C1_bed$end <- G6_C1_bed$start+1
G6_C1_bed$start <- G6_C1_bed$start*100000
G6_C1_bed$end <- G6_C1_bed$end*100000
G6_C1_bed$end <- as.numeric(G6_C1_bed$end)
G6_C1_bed$start <- as.numeric(G6_C1_bed$start)
G6_C1_bed <- select(G6_C1_bed,chr,start,end,meth_level)
G6_C1_bed <- na.omit(G6_C1_bed)
write.table(G6_C1_bed, "G6_C1.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_C2_bed <- as.data.frame(str_split_fixed(G6_C2$id,"_",2))
colnames(G6_C2_bed) <- c("chr","start")
G6_C2_bed$cite <- G6_C2$id
G6_C2_bed$meth_level <- G6_C2$`HepG2-SNHG6-C2`
G6_C2_bed$start <- as.numeric(G6_C2_bed$start)
G6_C2_bed$end <- G6_C2_bed$start+1
G6_C2_bed$start <- G6_C2_bed$start*100000
G6_C2_bed$end <- G6_C2_bed$end*100000
G6_C2_bed$end <- as.numeric(G6_C2_bed$end)
G6_C2_bed$start <- as.numeric(G6_C2_bed$start)
G6_C2_bed <- select(G6_C2_bed,chr,start,end,meth_level)
G6_C2_bed <- na.omit(G6_C2_bed)
write.table(G6_C2_bed, "G6_C2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_C3_bed <- as.data.frame(str_split_fixed(G6_C3$id,"_",2))
colnames(G6_C3_bed) <- c("chr","start")
G6_C3_bed$cite <- G6_C3$id
G6_C3_bed$meth_level <- G6_C3$`HepG2-SNHG6-C3`
G6_C3_bed$start <- as.numeric(G6_C3_bed$start)
G6_C3_bed$end <- G6_C3_bed$start+1
G6_C3_bed$start <- G6_C3_bed$start*100000
G6_C3_bed$end <- G6_C3_bed$end*100000
G6_C3_bed$end <- as.numeric(G6_C3_bed$end)
G6_C3_bed$start <- as.numeric(G6_C3_bed$start)
G6_C3_bed <- select(G6_C3_bed,chr,start,end,meth_level)
G6_C3_bed <- na.omit(G6_C3_bed)
write.table(G6_C3_bed, "G6_C3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_E1_bed <- as.data.frame(str_split_fixed(G6_E1$id,"_",2))
colnames(G6_E1_bed) <- c("chr","start")
G6_E1_bed$cite <- G6_E1$id
G6_E1_bed$meth_level <- G6_E1$`HepG2-SNHG6-OE1`
G6_E1_bed$start <- as.numeric(G6_E1_bed$start)
G6_E1_bed$end <- G6_E1_bed$start+1
G6_E1_bed$start <- G6_E1_bed$start*100000
G6_E1_bed$end <- G6_E1_bed$end*100000
G6_E1_bed$end <- as.numeric(G6_E1_bed$end)
G6_E1_bed$start <- as.numeric(G6_E1_bed$start)
G6_E1_bed <- select(G6_E1_bed,chr,start,end,meth_level)
G6_E1_bed <- na.omit(G6_E1_bed)
write.table(G6_E1_bed, "G6_E1.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_E2_bed <- as.data.frame(str_split_fixed(G6_E2$id,"_",2))
colnames(G6_E2_bed) <- c("chr","start")
G6_E2_bed$cite <- G6_E2$id
G6_E2_bed$meth_level <- G6_E2$`HepG2-SNHG6-OE2`
G6_E2_bed$start <- as.numeric(G6_E2_bed$start)
G6_E2_bed$end <- G6_E2_bed$start+1
G6_E2_bed$start <- G6_E2_bed$start*100000
G6_E2_bed$end <- G6_E2_bed$end*100000
G6_E2_bed$end <- as.numeric(G6_E2_bed$end)
G6_E2_bed$start <- as.numeric(G6_E2_bed$start)
G6_E2_bed <- select(G6_E2_bed,chr,start,end,meth_level)
G6_E2_bed <- na.omit(G6_E2_bed)
write.table(G6_E2_bed, "G6_E2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

G6_E3_bed <- as.data.frame(str_split_fixed(G6_E3$id,"_",2))
colnames(G6_E3_bed) <- c("chr","start")
G6_E3_bed$cite <- G6_E3$id
G6_E3_bed$meth_level <- G6_E3$`HepG2-SNHG6-OE3`
G6_E3_bed$start <- as.numeric(G6_E3_bed$start)
G6_E3_bed$end <- G6_E3_bed$start+1
G6_E3_bed$start <- G6_E3_bed$start*100000
G6_E3_bed$end <- G6_E3_bed$end*100000
G6_E3_bed$end <- as.numeric(G6_E3_bed$end)
G6_E3_bed$start <- as.numeric(G6_E3_bed$start)
G6_E3_bed <- select(G6_E3_bed,chr,start,end,meth_level)
G6_E3_bed <- na.omit(G6_E3_bed)
write.table(G6_E3_bed, "G6_E3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
