setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240705/bismark_cov/dmr_dmcp0.01delta0.2")
library(Rcircos)

install.packages("RCircos")
dmr_sm <-  fread("GA45_dmrs.bedGraph")
table(dmr_sm$V10)
#total
#Down-regulated   Up-regulated 
#79            117 

dmr_pmd <- fread("dmr_pmd.bedGraph")

dmr_pmd_sub <- subset(dmr_pmd,subset=V6!="Neither")

dmr_pmd_sub_p <- subset(dmr_pmd_sub,subset=V6=="commonPMD")

dmr_pmd_sub_h <- subset(dmr_pmd_sub,subset=V6=="commonHMD")

table(dmr_pmd_sub_h$V16)
#hmd
#Down-regulated   Up-regulated 
#14             12 

table(dmr_pmd_sub_p$V16)
#pmd
#Down-regulated   Up-regulated 
#24             42 

#pvh
sig_pmd_pvh_fisher_data <- matrix(c( 42,12,24,14), nrow = 2)
colnames(sig_pmd_pvh_fisher_data) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")
fisher.test(sig_pmd_pvh_fisher_data)
mcnemar.test(sig_pmd_pvh_fisher_data, correct=T) 
#p-value = 0.1599 odds ratio 2.02533
#McNemar's chi-squared = 3.3611, df = 1, p-value = 0.06675


#hmd
#Down-regulated   Up-regulated 
#79            117 
#Down-regulated   Up-regulated 
#14             12  
sig_hmd_fisher_data <- matrix(c(12,105,14,65), nrow = 2)
colnames(sig_hmd_fisher_data) <- c("up", "down") 
rownames(sig_hmd_fisher_data) <- c("hmd", "other")
fisher.test(sig_hmd_fisher_data)
chisq.test(sig_hmd_fisher_data,correct = F)
sig_hmd_chisq_result <- as.data.frame(unlist(chisq.test(sig_hmd_fisher_data)))
sig_hmd_fisher_result <- as.data.frame(unlist(fisher.test(sig_hmd_fisher_data)))
sig_hmd_result <- matrix(c(sig_hmd_chisq_result[1:3,],sig_hmd_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_hmd_result) <- 'hmd'
row.names(sig_hmd_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#odds ratio   0.53238 p-value = 0.1392
#X-squared = 2.2841, df = 1, p-value = 0.1307


#pmd
#Down-regulated   Up-regulated 
#79            117 
#Down-regulated   Up-regulated 
#24             42 
sig_pmd_fisher_data <- matrix(c(42,75,24,55), nrow = 2)
colnames(sig_pmd_fisher_data) <- c("up", "down") 
rownames(sig_pmd_fisher_data) <- c("pmd", "other")
fisher.test(sig_pmd_fisher_data)
chisq.test(sig_pmd_fisher_data,correct = F)
sig_pmd_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_fisher_data)))
sig_pmd_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_fisher_data)))
sig_pmd_result <- matrix(c(sig_pmd_chisq_result[1:3,],sig_pmd_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_result) <- 'pmd'
row.names(sig_pmd_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#odds ratio  1.281706 p-value = 0.4451
#X-squared = 0.64283, df = 1, p-value = 0.4227







dmr_H3K9me3 <- fread("dmr_H3K9me3.bedGraph")
dmr_H3K9me3_sub <- subset(dmr_H3K9me3,subset=V4>3)
dmr_H3K9me3_sub <- subset(dmr_sm,subset=V11 %in% unique(dmr_H3K9me3_sub$V15))
table(dmr_H3K9me3_sub$V10)
#total
#Down-regulated   Up-regulated 
#79            117 
#Down-regulated   Up-regulated 
#28             51 
sig_H3K9me3_fisher_data <- matrix(c(51,66,28,51), nrow = 2)
colnames(sig_H3K9me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K9me3_fisher_data) <- c("H3K9me3", "other")
fisher.test(sig_H3K9me3_fisher_data)
chisq.test(sig_H3K9me3_fisher_data,correct = F)
sig_H3K9me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_result <- matrix(c(sig_H3K9me3_chisq_result[1:3,],sig_H3K9me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K9me3_result) <- 'H3K9me3'
row.names(sig_H3K9me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#odds ratio  1.405006 p-value = 0.2994
#X-squared = 1.3008, df = 1, p-value = 0.2541











dmr_H3K27ac <- fread("dmr_H3K27ac.bedGraph")
dmr_H3K27ac_sub <- subset(dmr_H3K27ac,subset=V4>3)
dmr_H3K27ac_sub <- subset(dmr_sm,subset=V11 %in% unique(dmr_H3K27ac_sub$V15))
table(dmr_H3K27ac_sub$V10)
#total
#Down-regulated   Up-regulated 
#79            117 
#Down-regulated   Up-regulated 
# 36             48  
sig_H3K27ac_fisher_data <- matrix(c(48,69,36,43), nrow = 2)
colnames(sig_H3K27ac_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_fisher_data) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac_fisher_data)
chisq.test(sig_H3K27ac_fisher_data,correct = F)
sig_H3K27ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_result <- matrix(c(sig_H3K27ac_chisq_result[1:3,],sig_H3K27ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_result) <- 'H3K27ac'
row.names(sig_H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#odds ratio 0.831712  p-value = 0.5583
#X-squared = 0.3976, df = 1, p-value = 0.5283





dmr_H3K27ac_2 <- fread("dmr_H3K27ac_2.bedGraph")
table(dmr_H3K27ac_2$V20)
#total
#Down-regulated   Up-regulated 
#79            117 
#Down-regulated   Up-regulated 
#  19             25  
sig_H3K27ac_2_fisher_data <- matrix(c(25,92,19,60), nrow = 2)
colnames(sig_H3K27ac_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_2_fisher_data) <- c("H3K27ac_2", "other")
fisher.test(sig_H3K27ac_2_fisher_data)
chisq.test(sig_H3K27ac_2_fisher_data,correct = F)
sig_H3K27ac_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_result <- matrix(c(sig_H3K27ac_2_chisq_result[1:3,],sig_H3K27ac_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_2_result) <- 'H3K27ac_2'
row.names(sig_H3K27ac_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#odds ratio 0.8588034  p-value = 0.7279 
#X-squared = 0.19501, df = 1, p-value = 0.6588



dmr_H3K4me3 <- fread("dmr_H3K4me3.bedGraph")
dmr_H3K4me3_sub <- subset(dmr_H3K4me3,subset=V4>3)
dmr_H3K4me3_sub <- subset(dmr_sm,subset=V11 %in% unique(dmr_H3K4me3_sub$V15))
table(dmr_H3K4me3_sub$V10)
#total
#Down-regulated   Up-regulated 
#79            117 
#Down-regulated   Up-regulated 
# 46             57
sig_H3K4me3_fisher_data <- matrix(c(57,60,46,33), nrow = 2)
colnames(sig_H3K4me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_fisher_data) <- c("H3K4me3", "other")
fisher.test(sig_H3K4me3_fisher_data)
chisq.test(sig_H3K4me3_fisher_data,correct = F)
sig_H3K4me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_result <- matrix(c(sig_H3K4me3_chisq_result[1:3,],sig_H3K4me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_result) <- 'H3K4me3'
row.names(sig_H3K4me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#odds ratio  0.6828737  p-value = p-value = 0.2434
#X-squared = 1.7104, df = 1, p-value = 0.1909









dmr_H3K4me3_2 <- fread("dmr_H3K4me3_2.bedGraph")
dmr_H3K4me3_2_sub <- subset(dmr_H3K4me3_2,subset=V4>3)
dmr_H3K4me3_2_sub <- subset(dmr_sm,subset=V11 %in% unique(dmr_H3K4me3_2_sub$V15))
table(dmr_H3K4me3_2_sub$V10)
#total
#Down-regulated   Up-regulated 
#79            117 
#Down-regulated   Up-regulated 
#  26             29
sig_H3K4me3_2_fisher_data <- matrix(c(29,88,26,53), nrow = 2)
colnames(sig_H3K4me3_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_2_fisher_data) <- c("H3K4me3_2", "other")
fisher.test(sig_H3K4me3_2_fisher_data)
chisq.test(sig_H3K4me3_2_fisher_data,correct = F)
sig_H3K4me3_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_result <- matrix(c(sig_H3K4me3_2_chisq_result[1:3,],sig_H3K4me3_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_2_result) <- 'H3K4me3_2'
row.names(sig_H3K4me3_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#odds ratio  0.67318   p-value = 0.257
#X-squared = 1.5422, df = 1, p-value = 0.2143






dmr_sig_test_result <- as.data.frame(t(
  cbind(sig_hmd_result,sig_pmd_result,
        sig_H3K27ac_result,sig_H3K27ac_2_result,
        sig_H3K9me3_result,
        sig_H3K4me3_result,sig_H3K4me3_2_result)))
dmr_sig_test_result$histone_mark <- row.names(dmr_sig_test_result)

dmr_sig_test_result <- dplyr::select(dmr_sig_test_result,1,3,4,5,6)


dmr_sig_test_result$X_squared <- as.numeric(dmr_sig_test_result$X_squared)
dmr_sig_test_result$chisq_p <- as.numeric(dmr_sig_test_result$chisq_p)
dmr_sig_test_result$fisher_p <- as.numeric(dmr_sig_test_result$fisher_p)
dmr_sig_test_result$odds_ratio <- as.numeric(dmr_sig_test_result$odds_ratio)


dmr_sig_test_result[,"logp"] <- -log10(dmr_sig_test_result$chisq_p)
dmr_sig_test_result[,"lnR"] <- log(dmr_sig_test_result$odds_ratio)

dmr_sig_test_result <- dplyr::mutate(dmr_sig_test_result,case = case_when(lnR > 0 ~ 'Enriched in Hyper',
                                                                          lnR < 0 ~ 'Enriched in Hypo'))

dmr_sig_test_result$fdr <- p.adjust(dmr_sig_test_result$chisq_p, method  = "BH")
dmr_sig_test_result[,"logfdr"] <- -log10(dmr_sig_test_result$fdr)

ggplot(dmr_sig_test_result,aes(x=lnR,y=logfdr))+
  geom_point(size = 1,aes(color = case))+# 根据expression水平进行着色
  xlab(expression("ln odds ratio")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  theme_bw() +
  geom_text_repel(aes(label=histone_mark, color = case), size =5,hjust=1,vjust=1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Fisher test result of Histone Mark ")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
