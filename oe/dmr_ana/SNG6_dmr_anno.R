setwd("~/projects/hcc/analysis/oe/dss/SNG6_OE")
library(data.table)
library(ggrepel)
rm(list = ls())
save.image("SNG6_OE.Rdata")


dmr_sig <- fread("SNG6_dmrs.bedGraph")
dmr_sig$V11 <- paste(dmr_sig$V1,dmr_sig$V2,sep="_")

dmr_pmd <- fread("dmr_pmd.bedGraph")

dmr_pmd_sub <- subset(dmr_pmd,subset=V6!="Neither")

dmr_pmd_sub_p <- subset(dmr_pmd_sub,subset=V6=="commonPMD")

dmr_pmd_sub_h <- subset(dmr_pmd_sub,subset=V6=="commonHMD")


###total
table(dmr_sig$V10)
#Down-regulated   Up-regulated 
#   115             56

###pmd-hmd
table(dmr_pmd_sub_p$V16)
#Down-regulated   Up-regulated 
# 33             23  

table(dmr_pmd_sub_h$V16)
#Down-regulated   Up-regulated 
# 47             15


sig_pmd_p_fisher_data <- matrix(c(23, 33, 33, 82), nrow = 2)
colnames(sig_pmd_p_fisher_data) <- c("up", "down") 
rownames(sig_pmd_p_fisher_data) <- c("pmd", "other")
fisher.test(sig_pmd_p_fisher_data) 
chisq.test(sig_pmd_p_fisher_data)
sig_pmd_p_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_p_fisher_data)))
sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_p_fisher_data)))
sig_pmd_p_result <- matrix(c(sig_pmd_p_chisq_result[1:3,],sig_pmd_p_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_p_result) <- 'PMD'
row.names(sig_pmd_p_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.12 odds ratio 1.726137  
#X-squared = 2.0872, df = 1, p-value = 0.1485


sig_pmd_h_fisher_data <- matrix(c(15, 47, 41, 68), nrow = 2)
colnames(sig_pmd_h_fisher_data) <- c("up", "down") 
rownames(sig_pmd_h_fisher_data) <- c("hmd", "other")
fisher.test(sig_pmd_h_fisher_data) 
chisq.test(sig_pmd_h_fisher_data)
sig_pmd_h_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_h_fisher_data)))
sig_pmd_h_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_h_fisher_data)))
sig_pmd_h_result <- matrix(c(sig_pmd_h_chisq_result[1:3,],sig_pmd_h_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_h_result) <- 'HMD'
row.names(sig_pmd_h_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.09027  odds ratio   0.5312513
#X-squared = 2.6516, df = 1, p-value = 0.1034



sig_pmd_pvh_fisher_data <- matrix(c(23, 15, 33, 47), nrow = 2)
colnames(sig_pmd_pvh_fisher_data) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")

mcnemar.test(sig_pmd_pvh_fisher_data, correct=F) 
fisher.test(sig_pmd_pvh_fisher_data) 
#p-value = 0.07515  odds ratio 2.169114 
#McNemar's chi-squared = 6.75, df = 1, p-value = 0.009375












dmr_sig_H3K9me3 <- fread("dmr_H3K9me3.bedGraph")
table(dmr_sig_H3K9me3$V14)
dmr_sig_H3K9me3_sub <- subset(dmr_sig_H3K9me3,subset=V4>3)
dmr_sig_H3K9me3_sub$region <- paste(dmr_sig_H3K9me3_sub$V5,dmr_sig_H3K9me3_sub$V6,sep="_")
dmr_sig_H3K9me3_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K9me3_sub$region))
table(dmr_sig_H3K9me3_sub$V10)
#Down-regulated   Up-regulated 
#   115             56
#Down-regulated   Up-regulated 
#   35             17
sig_H3K9me3_fisher_data <- matrix(c(17, 39, 35, 80), nrow = 2)
colnames(sig_H3K9me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K9me3_fisher_data) <- c("H3K9me3", "other")
fisher.test(sig_H3K9me3_fisher_data) 
chisq.test(sig_H3K9me3_fisher_data)
sig_H3K9me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_result <- matrix(c(sig_H3K9me3_chisq_result[1:3,],sig_H3K9me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K9me3_result) <- 'H3K9me3'
row.names(sig_H3K9me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value =  1 odds ratio 0.9963583  
#X-squared = 0, df = 1, p-value = 1






dmr_sig_H3K27ac_2 <- fread("dmr_H3K27ac_2.bedGraph")
table(dmr_sig_H3K27ac_2$V20)
#Down-regulated   Up-regulated 
#    115             56
#Down-regulated   Up-regulated 
#    44             10 
sig_H3K27ac_2_fisher_data <- matrix(c(10, 46,44, 71), nrow = 2)
colnames(sig_H3K27ac_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_2_fisher_data) <- c("H3K27ac_2", "other")
fisher.test(sig_H3K27ac_2_fisher_data) 
chisq.test(sig_H3K27ac_2_fisher_data)
sig_H3K27ac_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_result <- matrix(c(sig_H3K27ac_2_chisq_result[1:3,],sig_H3K27ac_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_2_result) <- 'H3K27ac_2'
row.names(sig_H3K27ac_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = p-value = 0.008329 odds ratio  0.3528172 
#X-squared = 6.3428, df = 1, p-value = 0.01179


dmr_sig_H3K27ac <- fread("dmr_H3K27ac.bedGraph")
table(dmr_sig_H3K27ac$V14)
dmr_sig_H3K27ac_sub <- subset(dmr_sig_H3K27ac,subset=V4>3)
dmr_sig_H3K27ac_sub$region <- paste(dmr_sig_H3K27ac_sub$V5,dmr_sig_H3K27ac_sub$V6,sep="_")
dmr_sig_H3K27ac_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K27ac_sub$region))
table(dmr_sig_H3K27ac_sub$V10)
#Down-regulated   Up-regulated 
#    115             56
#Down-regulated   Up-regulated 
#   60             13
sig_H3K27ac_fisher_data <- matrix(c(13, 43,60, 55), nrow = 2)
colnames(sig_H3K27ac_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_fisher_data) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac_fisher_data) 
chisq.test(sig_H3K27ac_fisher_data)
sig_H3K27ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_result <- matrix(c(sig_H3K27ac_chisq_result[1:3,],sig_H3K27ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_result) <- 'H3K27ac'
row.names(sig_H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.0004811 odds ratio   0.2792337 
#X-squared = 11.753, df = 1, p-value = 0.0006074



dmr_sig_H3K4me3 <- fread("dmr_H3K4me3.bedGraph")
table(dmr_sig_H3K4me3$V14)
dmr_sig_H3K4me3_sub <- subset(dmr_sig_H3K4me3,subset=V4>3)
dmr_sig_H3K4me3_sub$region <- paste(dmr_sig_H3K4me3_sub$V5,dmr_sig_H3K4me3_sub$V6,sep="_")
dmr_sig_H3K4me3_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K4me3_sub$region))
table(dmr_sig_H3K4me3_sub$V10)
#Down-regulated   Up-regulated 
#   115             56 
#Down-regulated   Up-regulated 
#   69             20
sig_H3K4me3_fisher_data <- matrix(c(20, 36,69,46), nrow = 2)
colnames(sig_H3K4me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_fisher_data) <- c("H3K4me3", "other")
fisher.test(sig_H3K4me3_fisher_data) 
chisq.test(sig_H3K4me3_fisher_data)
sig_H3K4me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_result <- matrix(c(sig_H3K4me3_chisq_result[1:3,],sig_H3K4me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_result) <- 'H3K4me3'
row.names(sig_H3K4me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.003373 odds ratio   0.3725999  
#X-squared = 7.9533, df = 1, p-value = 0.0048


dmr_sig_H3K4me3_2 <- fread("dmr_H3K4me3_2.bedGraph")
table(dmr_sig_H3K4me3_2$V14)
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig_H3K4me3_2,subset=V4>3)
dmr_sig_H3K4me3_2_sub$region <- paste(dmr_sig_H3K4me3_2_sub$V5,dmr_sig_H3K4me3_2_sub$V6,sep="_")
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K4me3_2_sub$region))
table(dmr_sig_H3K4me3_2_sub$V10)
#Down-regulated   Up-regulated 
#      115             56 
#Down-regulated   Up-regulated 
#     47             11 
sig_H3K4me3_2_fisher_data <- matrix(c(11, 48,45, 68), nrow = 2)
colnames(sig_H3K4me3_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_2_fisher_data) <- c("H3K4me3_2", "other")
fisher.test(sig_H3K4me3_2_fisher_data) 
chisq.test(sig_H3K4me3_2_fisher_data)
sig_H3K4me3_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_result <- matrix(c(sig_H3K4me3_2_chisq_result[1:3,],sig_H3K4me3_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_2_result) <- 'H3K4me3_2'
row.names(sig_H3K4me3_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.00589 odds ratio  0.348336
#X-squared = 6.983, df = 1, p-value = 0.008229





###promoter
dmr_sig_promoter <- fread("dmr_promoter.bedGraph")
dmr_sig_promoter$region <- paste(dmr_sig_promoter$V7,dmr_sig_promoter$V8,sep="_")
dmr_sig_promoter_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_promoter$region))
table(dmr_sig_promoter_sub$V10)
#Down-regulated   Up-regulated 
# 47              9
#total
#Down-regulated   Up-regulated 
#   115             56 

dmr_sig_promoter_fisher_data <- matrix(c(9, 47,47,68), nrow = 2)
colnames(dmr_sig_promoter_fisher_data) <- c("up", "down") 
rownames(dmr_sig_promoter_fisher_data) <- c("promoter", "other")
fisher.test(dmr_sig_promoter_fisher_data)
chisq.test(dmr_sig_promoter_fisher_data)
sig_promoter_chisq_result <- as.data.frame(unlist(chisq.test(dmr_sig_promoter_fisher_data)))
sig_promoter_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_promoter_fisher_data)))
sig_promoter_result <- matrix(c(sig_promoter_chisq_result[1:3,],sig_promoter_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_promoter_result) <- 'promoter'
row.names(sig_promoter_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.001018 odds ratio 0.2789975
#X-squared = 9.4198, df = 1, p-value = 0.002147




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
  geom_text_repel(aes(label=histone_mark, color = case), size =5,hjust=1,vjust=-1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Chi-square test result of Histone Mark in SNHG6-OE-HepG2")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

