library(FactoMineR)
library(ggpubr)
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_202405/methy_mtx")
save.image("methy_mtx.Rdata")




hep_100k_e1_fragment <- hep_100k_e1$count[hep_100k_e1$count <= 5000]
head(hep_100k_e1_fragment)
hep_100k_e1_fragment <- data.frame(hep_100k_e1_fragment)
hep_100k_e1_res <- hist(hep_100k_e1_fragment$hep_100k_e1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, hep_100k_e1_res$breaks),
     y = c(0, 0, hep_100k_e1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "hep_100k_e1 Sample reads stat")

hep_100k_e2_fragment <- hep_100k_e2$count[hep_100k_e2$count <= 5000]
head(hep_100k_e2_fragment)
hep_100k_e2_fragment <- data.frame(hep_100k_e2_fragment)
hep_100k_e2_res <- hist(hep_100k_e2_fragment$hep_100k_e2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, hep_100k_e2_res$breaks),
     y = c(0, 0, hep_100k_e2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "hep_100k_e2 Sample reads stat")

hep_100k_e3_fragment <- hep_100k_e3$count[hep_100k_e3$count <= 5000]
head(hep_100k_e3_fragment)
hep_100k_e3_fragment <- data.frame(hep_100k_e3_fragment)
hep_100k_e3_res <- hist(hep_100k_e3_fragment$hep_100k_e3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, hep_100k_e3_res$breaks),
     y = c(0, 0, hep_100k_e3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "hep_100k_e3 Sample reads stat")

hep_100k_c1_fragment <- hep_100k_c1$count[hep_100k_c1$count <= 5000]
head(hep_100k_c1_fragment)
hep_100k_c1_fragment <- data.frame(hep_100k_c1_fragment)
hep_100k_c1_res <- hist(hep_100k_c1_fragment$hep_100k_c1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, hep_100k_c1_res$breaks),
     y = c(0, 0, hep_100k_c1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "hep_100k_c1 Sample reads stat")

hep_100k_c2_fragment <- hep_100k_c2$count[hep_100k_c2$count <= 5000]
head(hep_100k_c2_fragment)
hep_100k_c2_fragment <- data.frame(hep_100k_c2_fragment)
hep_100k_c2_res <- hist(hep_100k_c2_fragment$hep_100k_c2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, hep_100k_c2_res$breaks),
     y = c(0, 0, hep_100k_c2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "hep_100k_c2 Sample reads stat")

hep_100k_c3_fragment <- hep_100k_c3$count[hep_100k_c3$count <= 5000]
head(hep_100k_c3_fragment)
hep_100k_c3_fragment <- data.frame(hep_100k_c3_fragment)
hep_100k_c3_res <- hist(hep_100k_c3_fragment$hep_100k_c3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, hep_100k_c3_res$breaks),
     y = c(0, 0, hep_100k_c3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "hep_100k_c3 Sample reads stat")


hep_100k_c1_filt <- subset(hep_100k_c1,subset = hep_100k_c1$count >500)
hep_100k_c2_filt <- subset(hep_100k_c2,subset = hep_100k_c2$count >500)
hep_100k_c3_filt <- subset(hep_100k_c3,subset = hep_100k_c3$count >500)
hep_100k_e1_filt <- subset(hep_100k_e1,subset = hep_100k_e1$count >500)
hep_100k_e2_filt <- subset(hep_100k_e2,subset = hep_100k_e2$count >500)
hep_100k_e3_filt <- subset(hep_100k_e3,subset = hep_100k_e2$count >500)

hep_100k_c1_filt <- hep_100k_c1_filt[,c(2,6)]
hep_100k_c2_filt <- hep_100k_c2_filt[,c(2,6)]
hep_100k_c3_filt <- hep_100k_c3_filt[,c(2,6)]
hep_100k_e1_filt <- hep_100k_e1_filt[,c(2,6)]
hep_100k_e2_filt <- hep_100k_e2_filt[,c(2,6)]
hep_100k_e3_filt <- hep_100k_e3_filt[,c(2,6)]

hep_100k_id <- intersect(hep_100k_c1_filt$id,hep_100k_c2_filt$id)
hep_100k_id <- intersect(hep_100k_id,hep_100k_c3_filt$id)
hep_100k_id <- intersect(hep_100k_id,hep_100k_e1_filt$id)
hep_100k_id <- intersect(hep_100k_id,hep_100k_e2_filt$id)
hep_100k_id <- intersect(hep_100k_id,hep_100k_e3_filt$id)

hep_100k_c1_filt2 <- subset(hep_100k_c1_filt, subset= hep_100k_c1_filt$id %in% hep_100k_id)
hep_100k_c2_filt2 <- subset(hep_100k_c2_filt, subset= hep_100k_c2_filt$id %in% hep_100k_id)
hep_100k_c3_filt2 <- subset(hep_100k_c3_filt, subset= hep_100k_c3_filt$id %in% hep_100k_id)
hep_100k_e1_filt2 <- subset(hep_100k_e1_filt, subset= hep_100k_e1_filt$id %in% hep_100k_id)
hep_100k_e2_filt2 <- subset(hep_100k_e2_filt, subset= hep_100k_e2_filt$id %in% hep_100k_id)
hep_100k_e3_filt2 <- subset(hep_100k_e3_filt, subset= hep_100k_e3_filt$id %in% hep_100k_id)

hep_100k_filt <- cbind(hep_100k_c1_filt2,hep_100k_c2_filt2[,2],hep_100k_c3_filt2[,2],hep_100k_e1_filt2[,2],hep_100k_e2_filt2[,2],hep_100k_e3_filt2[,2])
hep_100k_filt <- as.data.frame(hep_100k_filt[,2:7])
row.names(hep_100k_filt) <- hep_100k_c1_filt2$id
hep_100k_filt <- na.omit(hep_100k_filt)
pheatmap(hep_100k_filt,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,annotation_row = pmd_row_anno,
         fontsize_col= 12,
         angle_col = 45,
         border=FALSE) 


hep_100k_baseMean = as.data.frame(colMeans(hep_100k_filt) )
hep_100k_PMDMean = as.data.frame(colMeans(na.omit(hep_100k_filt[pmd_regions,]))) 

hep_100k_baseMean$region <- "base"
hep_100k_PMDMean$region <- "pmd"
colnames(hep_100k_PMDMean) <- "mean"
colnames(hep_100k_baseMean) <- "mean"

hep_100k_baseMean$sample_type <- c(rep("control",3),rep("experiment",3))
hep_100k_PMDMean$sample_type <- c(rep("control",3),rep("experiment",3))


hep_100k_meansum <- rbind(hep_100k_PMDMean,hep_100k_baseMean)
hep_100k_meansum$sample_type <- c(rep("control",3),rep("experiment",3),rep("control",3),rep("experiment",3))
colnames(hep_100k_meansum) <- c("mean","region","sample_type")
ggplot(hep_100k_meansum,aes(x=sample_type,y=mean,fill=region))+
  geom_boxplot(width=0.5)+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("control","experiment")),
                     method = "t.test",label = "p.signif",
                     label.y =0.45 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(hep_100k_PMDMean,aes(x=sample_type,y=mean,fill=sample_type))+
    geom_boxplot(width=0.5)+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
    stat_compare_means(comparisons = list(c("control","experiment")),
                       method = "t.test",label = "p.signif",
                       label.y =0.24 )+theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),  ##去掉背景网格
          axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
          axis.title.x = element_blank(),
          legend.position = "none",
          legend.direction = "vertical",
          legend.title =element_blank())


ggplot(hep_100k_baseMean,aes(x=sample_type,y=mean,fill=sample_type))+
  geom_boxplot(width=0.5)+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("control","experiment")),
                     method = "t.test",label = "p.signif",
                     label.y =0.4 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hep_100k_t <- t(hep_100k_filt)
hep_100k_group <- data.frame(Sample = rownames(hep_10k_t), Group = rep(c("Control", "SH"), each = 3))
hep_100k_pca <- PCA(hep_100k_t, ncp = 2, scale.unit = TRUE, graph = FALSE)
hep_100k_pca_sample <- data.frame(hep_100k_pca$ind$coord[ ,1:2])
hep_100k_pca_sample$Sample=row.names(hep_100k_pca_sample)
hep_100k_pca_eig1 <- round(hep_10k_pca$eig[1,2], 2)
hep_100k_pca_eig2 <- round(hep_10k_pca$eig[2,2],2 )


hep_100k_pca_sample <- merge(hep_100k_pca_sample,hep_100k_group,by="Sample")
head(hep_100k_pca_sample)
hep_100k_p <- ggplot(data = hep_100k_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 2) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', hep_100k_pca_eig1, '%'), y = paste('PCA2:', hep_100k_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="HepG2_pca")

##t-test
hep_100k_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[1:3], row[4:6])$p.value
})

hep_100k_dt <- apply(hep_100k_filt, 1, function(row) {
  mean(row[1:3])-mean(row[4:6])
})

hep_100k_adjusted_p_values <- p.adjust(hep_100k_p_values, method  = "fdr")
hep_100k_adjusted_p_values_b <- p.adjust(hep_100k_p_values, method  = "bonferroni")
hep_100k_test <- cbind(hep_100k_p_values,hep_100k_adjusted_p_values,hep_100k_adjusted_p_values_b,hep_100k_dt)
hep_100k_test <- as.data.frame(hep_100k_test)
colnames(hep_100k_test) <- c("pvalue","bf_adj","fdr","fc")
hep_100k_test_significant <- subset(hep_100k_test,subset =  pvalue< 0.05)
hep_100k_test_significant_up <- subset(hep_100k_test_significant,subset = fc<0)
hep_100k_test_significant_down <- subset(hep_100k_test_significant,subset = fc>0)



hep_100k_r1_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[1:2], row[4:5])$p.value
})

hep_100k_r2_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[1:2], row[5:6])$p.value
})

hep_100k_r3_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[1:2], row[c(4,6)])$p.value
})

hep_100k_r4_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[2:3], row[4:5])$p.value
})

hep_100k_r5_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[2:3], row[5:6])$p.value
})

hep_100k_r6_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[2:3], row[c(4,6)])$p.value
})

hep_100k_r7_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[c(1,3)], row[4:5])$p.value
})

hep_100k_r8_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[c(1,3)], row[5:6])$p.value
})

hep_100k_r9_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[c(1,3)], row[c(4,6)])$p.value
})


hep_100k_r10_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,6)])$p.value
})
hep_100k_r11_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,4)])$p.value
})
hep_100k_r12_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[c(1,6)], row[c(2,5)])$p.value
})
hep_100k_r13_p_values <- apply(hep_100k_filt, 1, function(row) {
  t.test(row[c(3,6)], row[c(2,5)])$p.value
})



hep_100k_r1_test <- as.data.frame(cbind(hep_100k_r1_p_values,row.names(hep_100k_filt)))
hep_100k_r2_test <- as.data.frame(cbind(hep_100k_r2_p_values,row.names(hep_100k_filt)))
hep_100k_r3_test <- as.data.frame(cbind(hep_100k_r3_p_values,row.names(hep_100k_filt)))
hep_100k_r4_test <- as.data.frame(cbind(hep_100k_r4_p_values,row.names(hep_100k_filt)))
hep_100k_r5_test <- as.data.frame(cbind(hep_100k_r5_p_values,row.names(hep_100k_filt)))
hep_100k_r6_test <- as.data.frame(cbind(hep_100k_r6_p_values,row.names(hep_100k_filt)))
hep_100k_r7_test <- as.data.frame(cbind(hep_100k_r7_p_values,row.names(hep_100k_filt)))
hep_100k_r8_test <- as.data.frame(cbind(hep_100k_r8_p_values,row.names(hep_100k_filt)))
hep_100k_r9_test <- as.data.frame(cbind(hep_100k_r9_p_values,row.names(hep_100k_filt)))

hep_100k_r10_test <- as.data.frame(cbind(hep_100k_r10_p_values,row.names(hep_100k_filt)))
hep_100k_r11_test <- as.data.frame(cbind(hep_100k_r11_p_values,row.names(hep_100k_filt)))
hep_100k_r12_test <- as.data.frame(cbind(hep_100k_r12_p_values,row.names(hep_100k_filt)))
hep_100k_r13_test <- as.data.frame(cbind(hep_100k_r13_p_values,row.names(hep_100k_filt)))

hep_100k_r1_test_sig <- subset(hep_100k_r1_test,hep_100k_r1_test$hep_100k_r1_p_values < 0.05)
hep_100k_r2_test_sig <- subset(hep_100k_r2_test,hep_100k_r2_test$hep_100k_r2_p_values < 0.05)
hep_100k_r3_test_sig <- subset(hep_100k_r3_test,hep_100k_r3_test$hep_100k_r3_p_values < 0.05)
hep_100k_r4_test_sig <- subset(hep_100k_r4_test,hep_100k_r4_test$hep_100k_r4_p_values < 0.05)
hep_100k_r5_test_sig <- subset(hep_100k_r5_test,hep_100k_r5_test$hep_100k_r5_p_values < 0.05)
hep_100k_r6_test_sig <- subset(hep_100k_r6_test,hep_100k_r6_test$hep_100k_r6_p_values < 0.05)
hep_100k_r7_test_sig <- subset(hep_100k_r7_test,hep_100k_r7_test$hep_100k_r7_p_values < 0.05)
hep_100k_r8_test_sig <- subset(hep_100k_r8_test,hep_100k_r8_test$hep_100k_r8_p_values < 0.05)
hep_100k_r9_test_sig <- subset(hep_100k_r9_test,hep_100k_r9_test$hep_100k_r9_p_values < 0.05)
hep_100k_r10_test_sig <- subset(hep_100k_r10_test,hep_100k_r10_test$hep_100k_r10_p_values < 0.05)
hep_100k_r11_test_sig <- subset(hep_100k_r11_test,hep_100k_r11_test$hep_100k_r11_p_values < 0.05)
hep_100k_r12_test_sig <- subset(hep_100k_r12_test,hep_100k_r12_test$hep_100k_r12_p_values < 0.05)
hep_100k_r13_test_sig <- subset(hep_100k_r13_test,hep_100k_r13_test$hep_100k_r13_p_values < 0.05)


hep_100k_diff <- as.data.frame(cbind(c(651,504,507,703,558,538,600,620),
                               c("DA","DA","DA","DA","Control","Control","Control","Control")))
colnames(hep_100k_diff) <- c("DMR_num","group")
hep_100k_diff$DMR_num <- as.numeric(hep_100k_diff$DMR_num)
ggplot(hep_100k_diff,aes(x=group,y=DMR_num,fill=group))+
  geom_boxplot(width=0.5)+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("DA","Control")),
                     method = "t.test",label = "p.signif",
                     label.y =700 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())
  geom_point(aes(x=group,y=DMR_num))
