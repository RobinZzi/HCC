rm(list = ls())
library(pheatmap)
library(dplyr)
library(edgeR)
library(ggpubr)
library(FactoMineR)
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240705/methy_mtx")
save.image("20240705.Rdata")


setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/methy_mtx")
pmd_anno <- fread("PMD_coordinates_hg38.bed")
pmd_anno <- mutate(pmd_anno, region = paste(V1,V2/100000,sep = '_'))
pmd_row_anno <- dplyr::select(pmd_anno,region,V6)
pmd_row_anno <- as.data.frame(pmd_row_anno[,-1])
row.names(pmd_row_anno) <- pmd_anno$region
colnames(pmd_row_anno) <- "region_type"


setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240705/methy_mtx/100kbin")
ga45_100k_e1 <- fread("GA45_CRi1_100kbin.txt")
ga45_100k_e2 <- fread("GA45_CRi2_100kbin.txt")
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240701/methy_mtx/100kbin")
ga45_100k_c1 <- fread("GA45_V1_100kbin.txt")
ga45_100k_c2 <- fread("GA45_V2_100kbin.txt")


ga45_100k <- cbind(ga45_100k_e1,ga45_100k_e2[,6],ga45_100k_c1[,6],ga45_100k_c2[,6])
ga45_100k_sum <- as.data.frame(ga45_100k[,6:9])
row.names(ga45_100k_sum) <- ga45_100k$id
ga45_100k_sum <- na.omit(ga45_100k_sum)
colnames(ga45_100k_sum) <- c("GA45_CR1","GA45_CR2","GA45_V1","GA45_V2")
pheatmap(ga45_100k_sum,
         cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 10,
         angle_col = 45,
         border=FALSE,annotation_row = pmd_row_anno) 

ga45_100k_baseMean = as.data.frame(colMeans(ga45_100k_sum) ) 





pmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "commonPMD"))
hmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "commonHMD"))
ga45_100k_pmd <- ga45_100k_sum[pmd_regions,]
ga45_100k_hmd <- ga45_100k_sum[hmd_regions,]
ga45_100k_pmd <- na.omit(ga45_100k_pmd)
ga45_100k_hmd <- na.omit(ga45_100k_hmd)

ga45_100k_pmd_baseMean = as.data.frame(colMeans(ga45_100k_pmd) ) 
ga45_100k_hmd_baseMean = as.data.frame(colMeans(ga45_100k_hmd) ) 
ga45_100k_mean <- cbind(ga45_100k_baseMean,ga45_100k_pmd_baseMean,ga45_100k_hmd_baseMean)
colnames(ga45_100k_mean) <- c("global","PMD","HMD")

ga45_100k_PMDMean = as.data.frame(colMeans(ga45_100k_pmd) ) 
colnames(ga45_100k_PMDMean) <- 'mean_methylation_level'
ga45_100k_PMDMean$region <- 'PMD'
ga45_100k_PMDMean$group <- c("KD","KD","CT","CT")

ga45_100k_mean <- as.data.frame(rbind(ga45_100k_baseMean,ga45_100k_PMDMean))

colnames(ga45_100k_mean) <- c("global","PMD")
ga45_100k_mean[,"group"] <- c("KD","KD","CT","CT","KD","KD","CT","CT")
ga45_100k_mean$mean_methylation_level <- as.numeric(ga45_100k_mean$mean_methylation_level)
ggplot(ga45_100k_PMDMean,aes(group,mean_methylation_level,fill=group))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("KD","CT")),
                     method = "wilcox.test",label = "p.signif")+theme_bw()

setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240701/methy_mtx/10kbin")

ga45_10k_e1 <- fread("GA45_CRi1_10kbin.txt")
ga45_10k_e2 <- fread("GA45_CRi2_10kbin.txt")
ga45_10k_c1 <- fread("GA45_V1_10kbin.txt")
ga45_10k_c2 <- fread("GA45_V2_10kbin.txt")


ga45_10k <- cbind(ga45_10k_e1,ga45_10k_e2[,6],ga45_10k_c1[,6],ga45_10k_c2[,6])
ga45_10k_sum <- as.data.frame(ga45_10k[,6:9])
row.names(ga45_10k_sum) <- ga45_10k$id
ga45_10k_sum <- na.omit(ga45_10k_sum)
colnames(ga45_10k_sum) <- c("GA45_CR1","GA45_CR2","GA45_V1","GA45_V2")
pheatmap(ga45_10k_sum,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 10,
         angle_col = 45,
         border=FALSE,annotation_row = pmd_row_anno) 

ga45_10k_baseMean = as.data.frame(colMeans(ga45_10k_sum) ) 

pmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "commonPMD"))

ga45_10k_pmd <- ga45_10k_sum[pmd_regions,]
ga45_10k_pmd <- na.omit(ga45_10k_pmd)
ga45_10k_PMDMean = as.data.frame(colMeans(ga45_10k_pmd) ) 
colnames(ga45_10k_PMDMean) <- "293T-PMD"



ga45_10k_e1_fragment <- ga45_10k_e1$count[ga45_10k_e1$count <= 1000]
head(ga45_10k_e1_fragment)
ga45_10k_e1_fragment <- data.frame(ga45_10k_e1_fragment)
ga45_10k_e1_res <- hist(ga45_10k_e1_fragment$ga45_10k_e1_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, ga45_10k_e1_res$breaks),
     y = c(0, 0, ga45_10k_e1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "ga45_10k_e1_res Sample reads stat")


ga45_10k_e2_fragment <- ga45_10k_e2$count[ga45_10k_e2$count <= 1000]
head(ga45_10k_e2_fragment)
ga45_10k_e2_fragment <- data.frame(ga45_10k_e2_fragment)
ga45_10k_e2_res <- hist(ga45_10k_e2_fragment$ga45_10k_e2_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, ga45_10k_e2_res$breaks),
     y = c(0, 0, ga45_10k_e2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "ga45_10k_e2_res Sample reads stat")

ga45_10k_c1_fragment <- ga45_10k_c1$count[ga45_10k_c1$count <= 1000]
head(ga45_10k_c1_fragment)
ga45_10k_c1_fragment <- data.frame(ga45_10k_c1_fragment)
ga45_10k_c1_res <- hist(ga45_10k_c1_fragment$ga45_10k_c1_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, ga45_10k_c1_res$breaks),
     y = c(0, 0, ga45_10k_c1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "ga45_10k_c1_res Sample reads stat")

ga45_10k_c2_fragment <- ga45_10k_c2$count[ga45_10k_c2$count <= 1000]
head(ga45_10k_c2_fragment)
ga45_10k_c2_fragment <- data.frame(ga45_10k_c2_fragment)
ga45_10k_c2_res <- hist(ga45_10k_c2_fragment$ga45_10k_c2_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, ga45_10k_c2_res$breaks),
     y = c(0, 0, ga45_10k_c2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "ga45_10k_c2_res Sample reads stat")






ga45_10k_c1_filt <- subset(ga45_10k_c1,subset = ga45_10k_c1$count >5)
ga45_10k_c2_filt <- subset(ga45_10k_c2,subset = ga45_10k_c2$count >5)

ga45_10k_e1_filt <- subset(ga45_10k_e1,subset = ga45_10k_e1$count >5)
ga45_10k_e2_filt <- subset(ga45_10k_e2,subset = ga45_10k_e2$count >5)


ga45_10k_c1_filt <- ga45_10k_c1_filt[,c(2,6)]
ga45_10k_c2_filt <- ga45_10k_c2_filt[,c(2,6)]
ga45_10k_e1_filt <- ga45_10k_e1_filt[,c(2,6)]
ga45_10k_e2_filt <- ga45_10k_e2_filt[,c(2,6)]






ga45_10k_id <- Reduce(intersect,list(ga45_10k_c1_filt$id,ga45_10k_c2_filt$id,ga45_10k_e1_filt$id,ga45_10k_e2_filt$id))




ga45_10k_c1_filt2 <- subset(ga45_10k_c1_filt, subset= ga45_10k_c1_filt$id %in% ga45_10k_id)
ga45_10k_c2_filt2 <- subset(ga45_10k_c2_filt, subset= ga45_10k_c2_filt$id %in% ga45_10k_id)
ga45_10k_e1_filt2 <- subset(ga45_10k_e1_filt, subset= ga45_10k_e1_filt$id %in% ga45_10k_id)
ga45_10k_e2_filt2 <- subset(ga45_10k_e2_filt, subset= ga45_10k_e2_filt$id %in% ga45_10k_id)


ga45_10k_filt <- cbind(ga45_10k_c1_filt2,ga45_10k_c2_filt2[,2],ga45_10k_e1_filt2[,2],ga45_10k_e2_filt2[,2])
ga45_10k_filt <- as.data.frame(ga45_10k_filt[,2:5])
row.names(ga45_10k_filt) <- ga45_10k_c1_filt2$id
ga45_10k_filt <- na.omit(ga45_10k_filt)
pheatmap(ga45_10k_filt,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,annotation_row = pmd_row_anno,
         fontsize_col= 12,
         angle_col = 45,
         border=FALSE) 


ga45_10k_baseMean = as.data.frame(colMeans(ga45_10k_filt) )
ga45_10k_PMDMean = as.data.frame(colMeans(na.omit(ga45_10k_filt[pmd_regions,]))) 
ga45_10k_baseMean$region <- "base"
ga45_10k_PMDMean$region <- "pmd"
colnames(ga45_10k_PMDMean) <- "mean"
colnames(ga45_10k_baseMean) <- "mean"

ga45_10k_baseMean

ga45_10k_mean <-as.data.frame(cbind(ga45_10k_baseMean[,'mean'],ga45_10k_PMDMean[,'mean']))

colnames(ga45_10k_mean) <- c('global','PMD')
rownames(ga45_10k_mean) <- c("GA45_V1","GA45_V2","GA45_CR1","GA45_CR2")





ga45_10k_t <- t(ga45_10k_filt)
ga45_10k_group <- data.frame(Sample = rownames(ga45_10k_t), Group = rep(c("Control", "SH"), each = 2))
ga45_10k_pca <- PCA(ga45_10k_t, ncp = 2, scale.unit = TRUE, graph = FALSE)
ga45_10k_pca_sample <- data.frame(ga45_10k_pca$ind$coord[ ,1:2])
ga45_10k_pca_sample$Sample=row.names(ga45_10k_pca_sample)
ga45_10k_pca_eig1 <- round(ga45_10k_pca$eig[1,2], 2)
ga45_10k_pca_eig2 <- round(ga45_10k_pca$eig[2,2],2 )


ga45_10k_pca_sample <- merge(ga45_10k_pca_sample,ga45_10k_group,by="Sample")
head(ga45_10k_pca_sample)
ggplot(data = ga45_10k_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 2) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', ga45_10k_pca_eig1, '%'), y = paste('PCA2:', ga45_10k_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="ga457_10k_pca")










data <- as.data.frame(cbind(c(8.66,9.89,1.07,0.95),
                            c(0.627,0.619,0.633,0.632),
                            c(0.498,0.49,0.504,0.503),
                            c(0.794,0.785,0.8,0.797),
                            c("V","V","CR","CR")))
colnames(data) <- c("FPKM","Global_Methy","PMD_Methy","HMD_Methy","Group")                      
data$FPKM <- as.numeric(data$FPKM)
data$Global_Methy <- as.numeric(data$Global_Methy)
data$PMD_Methy <- as.numeric(data$PMD_Methy)
data$HMD_Methy <- as.numeric(data$HMD_Methy)

library(ggplot2)

library(ggpubr)
ggplot(data,aes(Group,HMD_Methy,color=Group))+
  geom_boxplot(width=0.5)+
  geom_jitter(width = 0.1,shape = 20,size=2)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("V","CR")),
                     method = "t.test",label = "p.signif",
                     label.y =0.82 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())






data2 <- as.data.frame(cbind(c(0.627,0.619,0.633,0.632,0.498,0.49,0.504,0.503,0.794,0.785,0.8,0.797),
                             c("Global_Mean","Global_Mean","Global_Mean","Global_Mean",
                               "PMD_Mean","PMD_Mean","PMD_Mean","PMD_Mean",
                               "HMD_Mean","HMD_Mean","HMD_Mean","HMD_Mean"),
                             c("V","V","CR","CR","V","V","CR","CR","V","V","CR","CR")))
colnames(data2) <- c("Methy_level","Region","Group")                      
data2$Methy_level <- as.numeric(data2$Methy_level)


ggplot(data2,aes(Group,Methy_level,color=Group))+
  geom_boxplot(width=0.5)+
  geom_jitter(width = 0.1,shape = 20,size=2)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("V","CR")),
                     method = "t.test",label = "p.signif",
                     label.y =0.55 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())+
  facet_wrap(~ Region)

