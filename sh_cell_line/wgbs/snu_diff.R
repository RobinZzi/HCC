library(FactoMineR)
library(ggpubr)


save.image("methy_mtx.Rdata")






snu_100k_e1_fragment <- snu_100k_e1$count[snu_100k_e1$count <= 5000]
head(snu_100k_e1_fragment)
snu_100k_e1_fragment <- data.frame(snu_100k_e1_fragment)
snu_100k_e1_res <- hist(snu_100k_e1_fragment$snu_100k_e1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, snu_100k_e1_res$breaks),
     y = c(0, 0, snu_100k_e1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "snu_100k_e1 Sample reads stat")

snu_100k_e2_fragment <- snu_100k_e2$count[snu_100k_e2$count <= 5000]
head(snu_100k_e2_fragment)
snu_100k_e2_fragment <- data.frame(snu_100k_e2_fragment)
snu_100k_e2_res <- hist(snu_100k_e2_fragment$snu_100k_e2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, snu_100k_e2_res$breaks),
     y = c(0, 0, snu_100k_e2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "snu_100k_e2 Sample reads stat")

snu_100k_e3_fragment <- snu_100k_e3$count[snu_100k_e3$count <= 5000]
head(snu_100k_e3_fragment)
snu_100k_e3_fragment <- data.frame(snu_100k_e3_fragment)
snu_100k_e3_res <- hist(snu_100k_e3_fragment$snu_100k_e3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, snu_100k_e3_res$breaks),
     y = c(0, 0, snu_100k_e3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "snu_100k_e3 Sample reads stat")

snu_100k_c1_fragment <- snu_100k_c1$count[snu_100k_c1$count <= 5000]
head(snu_100k_c1_fragment)
snu_100k_c1_fragment <- data.frame(snu_100k_c1_fragment)
snu_100k_c1_res <- hist(snu_100k_c1_fragment$snu_100k_c1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, snu_100k_c1_res$breaks),
     y = c(0, 0, snu_100k_c1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "snu_100k_c1 Sample reads stat")

snu_100k_c2_fragment <- snu_100k_c2$count[snu_100k_c2$count <= 5000]
head(snu_100k_c2_fragment)
snu_100k_c2_fragment <- data.frame(snu_100k_c2_fragment)
snu_100k_c2_res <- hist(snu_100k_c2_fragment$snu_100k_c2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, snu_100k_c2_res$breaks),
     y = c(0, 0, snu_100k_c2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "snu_100k_c2 Sample reads stat")

snu_100k_c3_fragment <- snu_100k_c3$count[snu_100k_c3$count <= 5000]
head(snu_100k_c3_fragment)
snu_100k_c3_fragment <- data.frame(snu_100k_c3_fragment)
snu_100k_c3_res <- hist(snu_100k_c3_fragment$snu_100k_c3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, snu_100k_c3_res$breaks),
     y = c(0, 0, snu_100k_c3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "snu_100k_c3 Sample reads stat")


snu_100k_c1_filt <- subset(snu_100k_c1,subset = snu_100k_c1$count >500)
snu_100k_c2_filt <- subset(snu_100k_c2,subset = snu_100k_c2$count >500)
snu_100k_c3_filt <- subset(snu_100k_c3,subset = snu_100k_c3$count >500)
snu_100k_e1_filt <- subset(snu_100k_e1,subset = snu_100k_e1$count >500)
snu_100k_e2_filt <- subset(snu_100k_e2,subset = snu_100k_e2$count >500)
snu_100k_e3_filt <- subset(snu_100k_e3,subset = snu_100k_e2$count >500)

snu_100k_c1_filt <- snu_100k_c1_filt[,c(2,6)]
snu_100k_c2_filt <- snu_100k_c2_filt[,c(2,6)]
snu_100k_c3_filt <- snu_100k_c3_filt[,c(2,6)]
snu_100k_e1_filt <- snu_100k_e1_filt[,c(2,6)]
snu_100k_e2_filt <- snu_100k_e2_filt[,c(2,6)]
snu_100k_e3_filt <- snu_100k_e3_filt[,c(2,6)]

snu_100k_id <- intersect(snu_100k_c1_filt$id,snu_100k_c2_filt$id)
snu_100k_id <- intersect(snu_100k_id,snu_100k_c3_filt$id)
snu_100k_id <- intersect(snu_100k_id,snu_100k_e1_filt$id)
snu_100k_id <- intersect(snu_100k_id,snu_100k_e2_filt$id)
snu_100k_id <- intersect(snu_100k_id,snu_100k_e3_filt$id)

snu_100k_c1_filt2 <- subset(snu_100k_c1_filt, subset= snu_100k_c1_filt$id %in% snu_100k_id)
snu_100k_c2_filt2 <- subset(snu_100k_c2_filt, subset= snu_100k_c2_filt$id %in% snu_100k_id)
snu_100k_c3_filt2 <- subset(snu_100k_c3_filt, subset= snu_100k_c3_filt$id %in% snu_100k_id)
snu_100k_e1_filt2 <- subset(snu_100k_e1_filt, subset= snu_100k_e1_filt$id %in% snu_100k_id)
snu_100k_e2_filt2 <- subset(snu_100k_e2_filt, subset= snu_100k_e2_filt$id %in% snu_100k_id)
snu_100k_e3_filt2 <- subset(snu_100k_e3_filt, subset= snu_100k_e3_filt$id %in% snu_100k_id)

snu_100k_filt <- cbind(snu_100k_c1_filt2,snu_100k_c2_filt2[,2],snu_100k_c3_filt2[,2],snu_100k_e1_filt2[,2],snu_100k_e2_filt2[,2],snu_100k_e3_filt2[,2])
snu_100k_filt <- as.data.frame(snu_100k_filt[,2:7])
row.names(snu_100k_filt) <- snu_100k_c1_filt2$id
snu_100k_filt <- na.omit(snu_100k_filt)
pheatmap(snu_100k_filt,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,annotation_row = pmd_row_anno,
         fontsize_col= 12,
         angle_col = 45,
         border=FALSE) 


snu_100k_baseMean = as.data.frame(colMeans(snu_100k_filt) )
snu_100k_PMDMean = as.data.frame(colMeans(na.omit(snu_100k_filt[pmd_regions,]))) 

snu_100k_baseMean$region <- "base"
snu_100k_PMDMean$region <- "pmd"
colnames(snu_100k_PMDMean) <- "mean"
colnames(snu_100k_baseMean) <- "mean"

snu_100k_baseMean$sample_type <- c(rep("control",3),rep("experiment",3))
snu_100k_PMDMean$sample_type <- c(rep("control",3),rep("experiment",3))


snu_100k_meansum <- rbind(snu_100k_PMDMean,snu_100k_baseMean)
snu_100k_meansum$sample_type <- c(rep("control",3),rep("experiment",3),rep("control",3),rep("experiment",3))
colnames(snu_100k_meansum) <- c("mean","region","sample_type")
ggplot(snu_100k_meansum,aes(x=sample_type,y=mean,fill=region))+
  geom_boxplot(width=0.5)+geom_jitter(width = 0.1,shape = 21, colour = "black")+
scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("control","experiment")),
                     method = "t.test",label = "p.signif",
                     label.y =0.65 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(snu_100k_PMDMean,aes(x=sample_type,y=mean,fill=sample_type))+
  geom_boxplot(width=0.5)+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("control","experiment")),
                     method = "t.test",label = "p.signif",
                     label.y =0.5 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


ggplot(snu_100k_baseMean,aes(x=sample_type,y=mean,fill=sample_type))+
  geom_boxplot(width=0.5)+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("control","experiment")),
                     method = "t.test",label = "p.signif",
                     label.y =0.615 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

snu_100k_t <- t(snu_100k_filt)
snu_100k_group <- data.frame(Sample = rownames(snu_10k_t), Group = rep(c("Control", "SH"), each = 3))
snu_100k_pca <- PCA(snu_100k_t, ncp = 2, scale.unit = TRUE, graph = FALSE)
snu_100k_pca_sample <- data.frame(snu_100k_pca$ind$coord[ ,1:2])
snu_100k_pca_sample$Sample=row.names(snu_100k_pca_sample)
snu_100k_pca_eig1 <- round(snu_10k_pca$eig[1,2], 2)
snu_100k_pca_eig2 <- round(snu_10k_pca$eig[2,2],2 )


snu_100k_pca_sample <- merge(snu_100k_pca_sample,snu_100k_group,by="Sample")
head(snu_100k_pca_sample)
snu_100k_p <- ggplot(data = snu_100k_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 2) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', snu_100k_pca_eig1, '%'), y = paste('PCA2:', snu_100k_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="SNU449_pca")

##t-test
snu_100k_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[1:3], row[4:6])$p.value
})

snu_100k_dt <- apply(snu_100k_filt, 1, function(row) {
  mean(row[1:3])-mean(row[4:6])
})

snu_100k_adjusted_p_values <- p.adjust(snu_100k_p_values, method  = "fdr")
snu_100k_adjusted_p_values_b <- p.adjust(snu_100k_p_values, method  = "bonferroni")
snu_100k_test <- cbind(snu_100k_p_values,snu_100k_adjusted_p_values,snu_100k_adjusted_p_values_b,snu_100k_dt)
snu_100k_test <- as.data.frame(snu_100k_test)
colnames(snu_100k_test) <- c("pvalue","bf_adj","fdr","fc")
snu_100k_test_significant <- subset(snu_100k_test,subset =  pvalue< 0.05)
snu_100k_test_significant_up <- subset(snu_100k_test_significant,subset = fc<0)
snu_100k_test_significant_down <- subset(snu_100k_test_significant,subset = fc>0)




snu_100k_r1_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[1:2], row[4:5])$p.value
})

snu_100k_r2_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[1:2], row[5:6])$p.value
})

snu_100k_r3_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[1:2], row[c(4,6)])$p.value
})

snu_100k_r4_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[2:3], row[4:5])$p.value
})

snu_100k_r5_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[2:3], row[5:6])$p.value
})

snu_100k_r6_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[2:3], row[c(4,6)])$p.value
})

snu_100k_r7_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[c(1,3)], row[4:5])$p.value
})

snu_100k_r8_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[c(1,3)], row[5:6])$p.value
})

snu_100k_r9_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[c(1,3)], row[c(4,6)])$p.value
})


snu_100k_r10_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,6)])$p.value
})
snu_100k_r11_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,4)])$p.value
})
snu_100k_r12_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[c(1,6)], row[c(2,5)])$p.value
})
snu_100k_r13_p_values <- apply(snu_100k_filt, 1, function(row) {
  t.test(row[c(3,6)], row[c(2,5)])$p.value
})



snu_100k_r1_test <- as.data.frame(cbind(snu_100k_r1_p_values,row.names(snu_100k_filt)))
snu_100k_r2_test <- as.data.frame(cbind(snu_100k_r2_p_values,row.names(snu_100k_filt)))
snu_100k_r3_test <- as.data.frame(cbind(snu_100k_r3_p_values,row.names(snu_100k_filt)))
snu_100k_r4_test <- as.data.frame(cbind(snu_100k_r4_p_values,row.names(snu_100k_filt)))
snu_100k_r5_test <- as.data.frame(cbind(snu_100k_r5_p_values,row.names(snu_100k_filt)))
snu_100k_r6_test <- as.data.frame(cbind(snu_100k_r6_p_values,row.names(snu_100k_filt)))
snu_100k_r7_test <- as.data.frame(cbind(snu_100k_r7_p_values,row.names(snu_100k_filt)))
snu_100k_r8_test <- as.data.frame(cbind(snu_100k_r8_p_values,row.names(snu_100k_filt)))
snu_100k_r9_test <- as.data.frame(cbind(snu_100k_r9_p_values,row.names(snu_100k_filt)))

snu_100k_r10_test <- as.data.frame(cbind(snu_100k_r10_p_values,row.names(snu_100k_filt)))
snu_100k_r11_test <- as.data.frame(cbind(snu_100k_r11_p_values,row.names(snu_100k_filt)))
snu_100k_r12_test <- as.data.frame(cbind(snu_100k_r12_p_values,row.names(snu_100k_filt)))
snu_100k_r13_test <- as.data.frame(cbind(snu_100k_r13_p_values,row.names(snu_100k_filt)))

snu_100k_r1_test_sig <- subset(snu_100k_r1_test,snu_100k_r1_test$snu_100k_r1_p_values < 0.05)
snu_100k_r2_test_sig <- subset(snu_100k_r2_test,snu_100k_r2_test$snu_100k_r2_p_values < 0.05)
snu_100k_r3_test_sig <- subset(snu_100k_r3_test,snu_100k_r3_test$snu_100k_r3_p_values < 0.05)
snu_100k_r4_test_sig <- subset(snu_100k_r4_test,snu_100k_r4_test$snu_100k_r4_p_values < 0.05)
snu_100k_r5_test_sig <- subset(snu_100k_r5_test,snu_100k_r5_test$snu_100k_r5_p_values < 0.05)
snu_100k_r6_test_sig <- subset(snu_100k_r6_test,snu_100k_r6_test$snu_100k_r6_p_values < 0.05)
snu_100k_r7_test_sig <- subset(snu_100k_r7_test,snu_100k_r7_test$snu_100k_r7_p_values < 0.05)
snu_100k_r8_test_sig <- subset(snu_100k_r8_test,snu_100k_r8_test$snu_100k_r8_p_values < 0.05)
snu_100k_r9_test_sig <- subset(snu_100k_r9_test,snu_100k_r9_test$snu_100k_r9_p_values < 0.05)
snu_100k_r10_test_sig <- subset(snu_100k_r10_test,snu_100k_r10_test$snu_100k_r10_p_values < 0.05)
snu_100k_r11_test_sig <- subset(snu_100k_r11_test,snu_100k_r11_test$snu_100k_r11_p_values < 0.05)
snu_100k_r12_test_sig <- subset(snu_100k_r12_test,snu_100k_r12_test$snu_100k_r12_p_values < 0.05)
snu_100k_r13_test_sig <- subset(snu_100k_r13_test,snu_100k_r13_test$snu_100k_r13_p_values < 0.05)


snu_100k_diff <- as.data.frame(cbind(c(591,688,625,568,626,602,613,609),
                                     c("DA","DA","DA","DA","Control","Control","Control","Control")))
colnames(snu_100k_diff) <- c("DMR_num","group")
snu_100k_diff$DMR_num <- as.numeric(snu_100k_diff$DMR_num)
ggplot(snu_100k_diff,aes(x=group,y=DMR_num,fill=group))+
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
