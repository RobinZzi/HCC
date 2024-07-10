setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/dss/HepG2/dmr_anno")
library(data.table)
library(ggrepel)
rm(list = ls())
save.image("GA45_KD.Rdata")


dmr_sig <- fread("GA45_dmrs.bedGraph")


dmr_pmd <- fread("dmr_pmd.bedGraph")

dmr_pmd_sub <- subset(dmr_pmd,subset=V6!="Neither")

dmr_pmd_sub_p <- subset(dmr_pmd_sub,subset=V6=="commonPMD")

dmr_pmd_sub_h <- subset(dmr_pmd_sub,subset=V6=="commonHMD")


###total
table(dmr_sig$V10)
#Down-regulated   Up-regulated 
# 221            218  

###pmd-hmd
table(dmr_pmd_sub_p$V16)
#Down-regulated   Up-regulated 
# 47             56 

table(dmr_pmd_sub_h$V16)
#Down-regulated   Up-regulated 
#    94             87


sig_pmd_p_fisher_data <- matrix(c(56, 162, 47, 174), nrow = 2)
colnames(sig_pmd_p_fisher_data) <- c("up", "down") 
rownames(sig_pmd_p_fisher_data) <- c("pmd", "other")
fisher.test(sig_pmd_p_fisher_data) 
chisq.test(sig_pmd_p_fisher_data)
sig_pmd_p_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_p_fisher_data)))
sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_p_fisher_data)))
sig_pmd_p_result <- matrix(c(sig_pmd_p_chisq_result[1:3,],sig_pmd_p_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_p_result) <- 'PMD'
row.names(sig_pmd_p_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.311 odds ratio 1.27902  
#X-squared = 0.96102, df = 1, p-value = 0.3269


sig_pmd_h_fisher_data <- matrix(c(87,131 , 94, 127), nrow = 2)
colnames(sig_pmd_h_fisher_data) <- c("up", "down") 
rownames(sig_pmd_h_fisher_data) <- c("hmd", "other")
fisher.test(sig_pmd_h_fisher_data) 
chisq.test(sig_pmd_h_fisher_data)
sig_pmd_h_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_h_fisher_data)))
sig_pmd_h_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_h_fisher_data)))
sig_pmd_h_result <- matrix(c(sig_pmd_h_chisq_result[1:3,],sig_pmd_h_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_h_result) <- 'HMD'
row.names(sig_pmd_h_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.628  odds ratio  0.8974942
#X-squared = 0.21329, df = 1, p-value = 0.6442



sig_pmd_pvh_fisher_data <- matrix(c(56, 87, 47, 94), nrow = 2)
colnames(sig_pmd_pvh_fisher_data) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")
fisher.test(sig_pmd_pvh_fisher_data) 
mcnemar.test(sig_pmd_pvh_fisher_data, correct=F) 

#p-value = 0.3254  odds ratio  1.286198  
#McNemar's chi-squared = 11.94, df = 1, p-value = 0.0005493












dmr_sig_H3K9me3 <- fread("dmr_H3K9me3.bedGraph")
table(dmr_sig_H3K9me3$V14)
dmr_sig_H3K9me3_sub <- subset(dmr_sig_H3K9me3,subset=V4>3)
dmr_sig_H3K9me3_sub$region <- paste(dmr_sig_H3K9me3_sub$V5,dmr_sig_H3K9me3_sub$V6,sep="_")
dmr_sig_H3K9me3_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K9me3_sub$region))
table(dmr_sig_H3K9me3_sub$V10)
#Down-regulated   Up-regulated 
# 221            218  
#Down-regulated   Up-regulated 
# 45             52 
sig_H3K9me3_fisher_data <- matrix(c(52, 166, 45, 176), nrow = 2)
colnames(sig_H3K9me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K9me3_fisher_data) <- c("H3K9me3", "other")
fisher.test(sig_H3K9me3_fisher_data) 
chisq.test(sig_H3K9me3_fisher_data)
sig_H3K9me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_result <- matrix(c(sig_H3K9me3_chisq_result[1:3,],sig_H3K9me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K9me3_result) <- 'H3K9me3'
row.names(sig_H3K9me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value =  0.4211 odds ratio 1.224596  
#X-squared = 0.5875, df = 1, p-value = 0.4434






dmr_sig_H3K27ac_2 <- fread("dmr_H3K27ac_2.bedGraph")
table(dmr_sig_H3K27ac_2$V20)
#Down-regulated   Up-regulated 
#   221            218  
#Down-regulated   Up-regulated 
#      76             87
sig_H3K27ac_2_fisher_data <- matrix(c(87, 131,76, 145), nrow = 2)
colnames(sig_H3K27ac_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_2_fisher_data) <- c("H3K27ac_2", "other")
fisher.test(sig_H3K27ac_2_fisher_data) 
chisq.test(sig_H3K27ac_2_fisher_data)
sig_H3K27ac_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_result <- matrix(c(sig_H3K27ac_2_chisq_result[1:3,],sig_H3K27ac_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_2_result) <- 'H3K27ac_2'
row.names(sig_H3K27ac_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.2376 odds ratio  1.26638 
#X-squared = 1.2054, df = 1, p-value = 0.2723


dmr_sig_H3K27ac <- fread("dmr_H3K27ac.bedGraph")
table(dmr_sig_H3K27ac$V14)
dmr_sig_H3K27ac_sub <- subset(dmr_sig_H3K27ac,subset=V4>3)
dmr_sig_H3K27ac_sub$region <- paste(dmr_sig_H3K27ac_sub$V5,dmr_sig_H3K27ac_sub$V6,sep="_")
dmr_sig_H3K27ac_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K27ac_sub$region))
table(dmr_sig_H3K27ac_sub$V10)
#Down-regulated   Up-regulated 
#   221            218  
#Down-regulated   Up-regulated 
#     97            101 
sig_H3K27ac_fisher_data <- matrix(c(101, 117,97, 124), nrow = 2)
colnames(sig_H3K27ac_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_fisher_data) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac_fisher_data) 
chisq.test(sig_H3K27ac_fisher_data)
sig_H3K27ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_result <- matrix(c(sig_H3K27ac_chisq_result[1:3,],sig_H3K27ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_result) <- 'H3K27ac'
row.names(sig_H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.6322 odds ratio   1.103284 
#X-squared = 0.17434, df = 1, p-value = 0.6763



dmr_sig_H3K4me3 <- fread("dmr_H3K4me3.bedGraph")
table(dmr_sig_H3K4me3$V14)
dmr_sig_H3K4me3_sub <- subset(dmr_sig_H3K4me3,subset=V4>3)
dmr_sig_H3K4me3_sub$region <- paste(dmr_sig_H3K4me3_sub$V5,dmr_sig_H3K4me3_sub$V6,sep="_")
dmr_sig_H3K4me3_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K4me3_sub$region))
table(dmr_sig_H3K4me3_sub$V10)
#Down-regulated   Up-regulated 
#   221            218  
#Down-regulated   Up-regulated 
#  115            115 
sig_H3K4me3_fisher_data <- matrix(c(115, 103,115,106), nrow = 2)
colnames(sig_H3K4me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_fisher_data) <- c("H3K4me3", "other")
fisher.test(sig_H3K4me3_fisher_data) 
chisq.test(sig_H3K4me3_fisher_data)
sig_H3K4me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_result <- matrix(c(sig_H3K4me3_chisq_result[1:3,],sig_H3K4me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_result) <- 'H3K4me3'
row.names(sig_H3K4me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.9239 odds ratio  1.029071 
#X-squared = 0.0029856, df = 1, p-value = 0.9564


dmr_sig_H3K4me3_2 <- fread("dmr_H3K4me3_2.bedGraph")
table(dmr_sig_H3K4me3_2$V14)
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig_H3K4me3_2,subset=V4>3)
dmr_sig_H3K4me3_2_sub$region <- paste(dmr_sig_H3K4me3_2_sub$V5,dmr_sig_H3K4me3_2_sub$V6,sep="_")
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K4me3_2_sub$region))
table(dmr_sig_H3K4me3_2_sub$V10)
#Down-regulated   Up-regulated 
#    221            218 
#Down-regulated   Up-regulated 
#     81             90
sig_H3K4me3_2_fisher_data <- matrix(c(90, 128,81, 140), nrow = 2)
colnames(sig_H3K4me3_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_2_fisher_data) <- c("H3K4me3_2", "other")
fisher.test(sig_H3K4me3_2_fisher_data) 
chisq.test(sig_H3K4me3_2_fisher_data)
sig_H3K4me3_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_result <- matrix(c(sig_H3K4me3_2_chisq_result[1:3,],sig_H3K4me3_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_2_result) <- 'H3K4me3_2'
row.names(sig_H3K4me3_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.3293 odds ratio  1.214731 
#X-squared = 0.8053, df = 1, p-value = 0.3695





###promoter
dmr_sig_promoter <- fread("dmr_promoter.bedGraph")
dmr_sig_promoter$region <- paste(dmr_sig_promoter$V7,dmr_sig_promoter$V8,sep="_")
dmr_sig_promoter_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_promoter$region))
table(dmr_sig_promoter_sub$V10)
#Down-regulated   Up-regulated 
#221            218 
#total
#Down-regulated   Up-regulated 
#    81             72 

dmr_sig_promoter_fisher_data <- matrix(c(72, 146,81,140), nrow = 2)
colnames(dmr_sig_promoter_fisher_data) <- c("up", "down") 
rownames(dmr_sig_promoter_fisher_data) <- c("promoter", "other")
fisher.test(dmr_sig_promoter_fisher_data)
chisq.test(dmr_sig_promoter_fisher_data)
sig_promoter_chisq_result <- as.data.frame(unlist(chisq.test(dmr_sig_promoter_fisher_data)))
sig_promoter_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_promoter_fisher_data)))
sig_promoter_result <- matrix(c(sig_promoter_chisq_result[1:3,],sig_promoter_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_promoter_result) <- 'promoter'
row.names(sig_promoter_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.4833 odds ratio 0.8526721
#X-squared = 0.48523, df = 1, p-value = 0.4861




dmr_sig_test_result <- as.data.frame(t(
  cbind(sig_promoter_result,
        sig_H3K27ac_2_result,sig_H3K27ac_result,
        sig_H3K9me3_result,
        sig_H3K4me3_result,sig_H3K4me3_2_result,
        sig_pmd_p_result,sig_pmd_h_result)))
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
  geom_text_repel(aes(label=histone_mark, color = case), size =5,hjust=0.5,vjust=-1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Chi-square test result of Histone Mark in GADD45A-KD-HepG2")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

