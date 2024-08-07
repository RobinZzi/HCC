rm(list = ls())
library(pheatmap)
library(dplyr)
library(edgeR)
library(FactoMineR)
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240611/methy_mtx/100kbin")

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/methy_mtx")
pmd_anno <- fread("ga_global_PMD.bedGraph")
pmd_row_anno <- dplyr::select(pmd_anno,V4,V9)
pmd_row_anno <- as.data.frame(pmd_row_anno[,-1])
row.names(pmd_row_anno) <- pmd_anno$V4
colnames(pmd_row_anno) <- "region_type"
save.image("20240611.Rdata")


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

pmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "PMD"))

huh_100k_pmd <- huh_100k_sum[pmd_regions,]
huh_100k_pmd <- na.omit(huh_100k_pmd)
huh_100k_PMDMean = as.data.frame(colMeans(huh_100k_pmd) ) 
colnames(huh_100k_PMDMean) <- "Huh7-PMD"
hmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "HMD"))

huh_100k_hmd <- huh_100k_sum[hmd_regions,]
huh_100k_hmd <- na.omit(huh_100k_hmd)
huh_100k_HMDMean = as.data.frame(colMeans(huh_100k_hmd) ) 
colnames(huh_100k_HMDMean) <- "Huh7-HMD"
colnames(huh_100k_baseMean) <- "Huh7-global"

huh_meansum <- cbind(huh_100k_PMDMean,huh_100k_HMDMean,huh_100k_baseMean)


group = c("sh","sh","v","v")
design = model.matrix(~0+group)


y = DGEList(counts=huh_10k_sum)
y<-calcNormFactors(y)
y<-estimateCommonDisp(y, rowsum.filter=0)
y<-estimateGLMTagwiseDisp(y,design) #
pdf("MDSplot.pdf")
plotMDS(y,label=group,xlim=c(-8,8))
dev.off()


fit_tag<-glmFit(y,design)
lrt.tagwise<-glmLRT(fit_tag,contrast=c(1,-1))
pvals_tag <- lrt.tagwise$table$PValue
FDR_tag<- p.adjust(pvals_tag, method="BH")
out2 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
out2 = out2[order(out1$FDR_tag),]

table(out1$FDR_tag<0.05)
sig1 = out1[which(out1$FDR_tag<0.05),]



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

huh_10k_baseMean = as.data.frame(colMeans(huh_10k_sum) ) 


huh_10k_pmd <- huh_10k_sum[pmd_regions,]
huh_10k_pmd <- na.omit(huh_10k_pmd)
huh_10k_PMDMean = as.data.frame(colMeans(huh_10k_pmd) ) 
colnames(huh_10k_PMDMean) <- "Huh7-PMD"



huh_10k_e1_fragment <- huh_10k_e1$count[huh_10k_e1$count <= 100]
head(huh_10k_e1_fragment)
huh_10k_e1_fragment <- data.frame(huh_10k_e1_fragment)
huh_10k_e1_res <- hist(huh_10k_e1_fragment$huh_10k_e1_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, huh_10k_e1_res$breaks),
     y = c(0, 0, huh_10k_e1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "huh_10k_e1_res Sample reads stat")


huh_10k_e2_fragment <- huh_10k_e2$count[huh_10k_e2$count <= 100]
head(huh_10k_e2_fragment)
huh_10k_e2_fragment <- data.frame(huh_10k_e2_fragment)
huh_10k_e2_res <- hist(huh_10k_e2_fragment$huh_10k_e2_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, huh_10k_e2_res$breaks),
     y = c(0, 0, huh_10k_e2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "huh_10k_e2_res Sample reads stat")

huh_10k_c1_fragment <- huh_10k_c1$count[huh_10k_c1$count <= 100]
head(huh_10k_c1_fragment)
huh_10k_c1_fragment <- data.frame(huh_10k_c1_fragment)
huh_10k_c1_res <- hist(huh_10k_c1_fragment$huh_10k_c1_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, huh_10k_c1_res$breaks),
     y = c(0, 0, huh_10k_c1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "huh_10k_c1_res Sample reads stat")

huh_10k_c2_fragment <- huh_10k_c2$count[huh_10k_c2$count <= 100]
head(huh_10k_c2_fragment)
huh_10k_c2_fragment <- data.frame(huh_10k_c2_fragment)
huh_10k_c2_res <- hist(huh_10k_c2_fragment$huh_10k_c2_fragment, breaks = 100, plot = FALSE)
plot(x = c(0, huh_10k_c2_res$breaks),
     y = c(0, 0, huh_10k_c2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "huh_10k_c2_res Sample reads stat")






huh_10k_c1_filt <- subset(huh_10k_c1,subset = huh_10k_c1$count >5)
huh_10k_c2_filt <- subset(huh_10k_c2,subset = huh_10k_c2$count >5)

huh_10k_e1_filt <- subset(huh_10k_e1,subset = huh_10k_e1$count >5)
huh_10k_e2_filt <- subset(huh_10k_e2,subset = huh_10k_e2$count >5)


huh_10k_c1_filt <- huh_10k_c1_filt[,c(2,6)]
huh_10k_c2_filt <- huh_10k_c2_filt[,c(2,6)]
huh_10k_e1_filt <- huh_10k_e1_filt[,c(2,6)]
huh_10k_e2_filt <- huh_10k_e2_filt[,c(2,6)]






huh_10k_id <- Reduce(intersect,list(huh_10k_c1_filt$id,huh_10k_c2_filt$id,huh_10k_e1_filt$id,huh_10k_e2_filt$id))




huh_10k_c1_filt2 <- subset(huh_10k_c1_filt, subset= huh_10k_c1_filt$id %in% huh_10k_id)
huh_10k_c2_filt2 <- subset(huh_10k_c2_filt, subset= huh_10k_c2_filt$id %in% huh_10k_id)
huh_10k_e1_filt2 <- subset(huh_10k_e1_filt, subset= huh_10k_e1_filt$id %in% huh_10k_id)
huh_10k_e2_filt2 <- subset(huh_10k_e2_filt, subset= huh_10k_e2_filt$id %in% huh_10k_id)


huh_10k_filt <- cbind(huh_10k_c1_filt2,huh_10k_c2_filt2[,2],huh_10k_e1_filt2[,2],huh_10k_e2_filt2[,2])
huh_10k_filt <- as.data.frame(huh_10k_filt[,2:5])
row.names(huh_10k_filt) <- huh_10k_c1_filt2$id
huh_10k_filt <- na.omit(huh_10k_filt)
pheatmap(huh_10k_filt,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,annotation_row = pmd_row_anno,
         fontsize_col= 12,
         angle_col = 45,
         border=FALSE) 


huh_10k_baseMean = as.data.frame(colMeans(huh_10k_filt) )
huh_10k_PMDMean = as.data.frame(colMeans(na.omit(huh_10k_filt[pmd_regions,]))) 
huh_10k_baseMean$region <- "base"
huh_10k_PMDMean$region <- "pmd"
colnames(huh_10k_PMDMean) <- "mean"
colnames(huh_10k_baseMean) <- "mean"

huh_10k_p_values <- apply(huh_10k_filt, 1, function(row) {
  t.test(row[1:2], row[3:4])$p.value
})

huh_10k_dt <- apply(huh_10k_filt, 1, function(row) {
  mean(row[1:2])-mean(row[3:4])
})










huh_10k_t <- t(huh_10k_filt)
huh_10k_group <- data.frame(Sample = rownames(huh_10k_t), Group = rep(c("Control", "SH"), each = 2))
huh_10k_pca <- PCA(huh_10k_t, ncp = 2, scale.unit = TRUE, graph = FALSE)
huh_10k_pca_sample <- data.frame(huh_10k_pca$ind$coord[ ,1:2])
huh_10k_pca_sample$Sample=row.names(huh_10k_pca_sample)
huh_10k_pca_eig1 <- round(huh_10k_pca$eig[1,2], 2)
huh_10k_pca_eig2 <- round(huh_10k_pca$eig[2,2],2 )


huh_10k_pca_sample <- merge(huh_10k_pca_sample,huh_10k_group,by="Sample")
head(huh_10k_pca_sample)
ggplot(data = huh_10k_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 2) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', huh_10k_pca_eig1, '%'), y = paste('PCA2:', huh_10k_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="Huh7_10k_pca")




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


