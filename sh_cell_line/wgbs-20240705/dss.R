library(DSS)
library(ggrepel)
library(dplyr)
library(stringr)
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240705/bismark_cov")
rm(list = ls())
load("dss.Rdata")
save.image("dss.Rdata")

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_cov")
GA45_c1 <- fread("GA45_V1.bismark.cov.gz")
GA45_c2 <- fread("GA45_V2.bismark.cov.gz")
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240705/bismark_cov")
GA45_e1 <- fread("GA45_CRi1.bismark.cov.gz")
GA45_e2 <- fread("GA45_CRi2.bismark.cov.gz")


GA45_c1$total <- GA45_c1$V5+GA45_c1$V6
GA45_c1_con <- dplyr::select(GA45_c1,1,2,7,5)
colnames(GA45_c1_con) <- c("chr","pos","N","X")

GA45_c2$total <- GA45_c2$V5+GA45_c2$V6
GA45_c2_con <- dplyr::select(GA45_c2,1,2,7,5)
colnames(GA45_c2_con) <- c("chr","pos","N","X")

GA45_e1$total <- GA45_e1$V5+GA45_e1$V6
GA45_e1_con <- dplyr::select(GA45_e1,1,2,7,5)
colnames(GA45_e1_con) <- c("chr","pos","N","X")

GA45_e2$total <- GA45_e2$V5+GA45_e2$V6
GA45_e2_con <- dplyr::select(GA45_e2,1,2,7,5)
colnames(GA45_e2_con) <- c("chr","pos","N","X")

BSobj <- makeBSseqData( list(GA45_c1_con, GA45_c2_con, GA45_e1_con, GA45_e2_con),
                        c("C1","C2", "E1", "E2") )
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("E1", "E2"))

dmlTest.sm <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("E1", "E2"), smoothing=TRUE)
dmlTest.sm_sig <- subset(dmlTest.sm,subset=fdr<0.05)
dmr_sig <- callDMR(dmlTest.sm_sig)
dmr_sig2 <- callDMR(dmlTest.sm_sig,p.threshold=1)

dmrs <- callDMR(dmlTest.sm_sig, p.threshold=0.01,delt=0.2,minlen = 50, minCG = 3, dis.merge = 100)
showOneDMR(dmrs['712',], BSobj)


dmrs <- mutate(dmrs,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                            diff.Methy>0 ~ "Down-regulated", # 下调
                                            TRUE ~ "Unchanged"))
table(dmrs$state)

dmr_sig <- mutate(dmr_sig,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                                       diff.Methy>0 ~ "Down-regulated", # 下调
                                                       TRUE ~ "Unchanged"))

dmr_sig$region <- paste(dmr_sig$chr,dmr_sig$start,sep = "_")
dmrs$region <- paste(dmrs$chr,dmrs$start,sep = "_")

dml_sig <- callDML(dmlTest.sm_sig)
dml_sig <- mutate(dml_sig,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                              diff>0 ~ "Down-regulated", # 下调
                                              TRUE ~ "Unchanged"))
table(dml_sig$state)
table(dmr_sig$state)
showOneDMR(dmrs['361941',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['2408',], BSobj)

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_cov")
dmr_sig_light <- fread("GA45_dmrs_sig.bedGraph")
length(intersect(dmr_sig_light$V11,dmr_sig$region))
intersect(dmr_sig_light$V9,dmr_sig$areaStat)


write.table(dmr_sig, "GA45_dmrs_sig.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
write.table(dmrs, "GA45_dmrs.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


write.table(dmrs.sm, "GA45_dmrs.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

df_color <- c("#1f76b6", "#ff7d0e")
dmrs.sm.prop <- as.data.frame(t(table(dmrs.sm$state)))
ggplot(dmrs.sm.prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))

dmls.sm.prop <- as.data.frame(t(table(dmls.sm$state)))
ggplot(dmls.sm.prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))


dmrs_sig_prop <- as.data.frame(t(table(dmlTest.sm_sig_dmr$state)))
ggplot(dmrs_sig_prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))

dmls_sig_prop <- as.data.frame(t(table(dmls_sig$state)))
ggplot(dmrs_sig_prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))


PMD <- fread("PMD_coordinates_hg38.bed")

dmr_pmd <- fread("dmr_sig_pmd.bedGraph")

dmr_pmd_sub <- subset(dmr_pmd,subset=V6!="Neither")

dmr_pmd_sub_p <- subset(dmr_pmd_sub,subset=V6=="commonPMD")

dmr_pmd_sub_h <- subset(dmr_pmd_sub,subset=V6=="commonHMD")

dmrs_pmd <- fread("dmrs_pmd.bedGraph")

dmrs_pmd_sub <- subset(dmrs_pmd,subset=V6!="Neither")

dmrs_pmd_sub_p <- subset(dmrs_pmd_sub,subset=V6=="commonPMD")

dmrs_pmd_sub_h <- subset(dmrs_pmd_sub,subset=V6=="commonHMD")

###total
table(dmr_sig$state)
#Down-regulated   Up-regulated 
#  793           1593 

###pmd-hmd
table(dmr_pmd_sub_p$V16)
#Down-regulated   Up-regulated 
#  217            454 

table(dmr_pmd_sub_h$V16)
#Down-regulated   Up-regulated 
# 170            407


sig_pmd_p_fisher_data <- matrix(c(454, 1139, 217, 576), nrow = 2)
colnames(sig_pmd_p_fisher_data) <- c("up", "down") 
rownames(sig_pmd_p_fisher_data) <- c("pmd", "other")
fisher.test(sig_pmd_p_fisher_data) 
chisq.test(sig_pmd_p_fisher_data)
sig_pmd_p_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_p_fisher_data)))
sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_p_fisher_data)))
sig_pmd_p_result <- matrix(c(sig_pmd_p_chisq_result[1:3,],sig_pmd_p_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_p_result) <- 'PMD'
row.names(sig_pmd_p_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.595 odds ratio 1.058011 
#X-squared = 0.28374, df = 1, p-value = 0.5943


sig_pmd_h_fisher_data <- matrix(c(407, 1186, 170, 663), nrow = 2)
colnames(sig_pmd_h_fisher_data) <- c("up", "down") 
rownames(sig_pmd_h_fisher_data) <- c("hmd", "other")
fisher.test(sig_pmd_h_fisher_data) 
chisq.test(sig_pmd_h_fisher_data)
sig_pmd_h_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_h_fisher_data)))
sig_pmd_h_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_h_fisher_data)))
sig_pmd_h_result <- matrix(c(sig_pmd_h_chisq_result[1:3,],sig_pmd_h_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_h_result) <- 'HMD'
row.names(sig_pmd_h_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.004905  odds ratio 1.338205 
#X-squared = 7.6943, df = 1, p-value = 0.005539



sig_pmd_pvh_fisher_data <- matrix(c(454, 217, 407, 170), nrow = 2)
colnames(sig_pmd_pvh_fisher_data) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")

mcnemar.test(sig_pmd_pvh_fisher_data, correct=F) 

sig_pmd_pvh_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_pvh_fisher_data) ))
#p-value = 0.2968  odds ratio 0.8739731 
#McNemar's chi-squared = 57.853, df = 1, p-value = 2.825e-14












dmr_sig_H3K9me3 <- fread("dmr_sig_H3K9me3.bedGraph")
table(dmr_sig_H3K9me3$V14)
dmr_sig_H3K9me3_sub <- subset(dmr_sig_H3K9me3,subset=V4>3)
dmr_sig_H3K9me3_sub <- subset(dmr_sig_H3K9me3,subset=V4>5)
dmr_sig_H3K9me3_sub$region <- paste(dmr_sig_H3K9me3_sub$V5,dmr_sig_H3K9me3_sub$V6,sep="_")
dmr_sig_H3K9me3_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K9me3_sub$region))
table(dmr_sig_H3K9me3_sub$state)
#Down-regulated   Up-regulated 
#  793           1593 
#Down-regulated   Up-regulated 
#  285            598 
sig_H3K9me3_fisher_data <- matrix(c(203, 1390, 96, 697), nrow = 2)
colnames(sig_H3K9me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K9me3_fisher_data) <- c("H3K9me3", "other")
fisher.test(sig_H3K9me3_fisher_data) 
chisq.test(sig_H3K9me3_fisher_data)
sig_H3K9me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K9me3_fisher_data)))
sig_H3K9me3_result <- matrix(c(sig_H3K9me3_chisq_result[1:3,],sig_H3K9me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K9me3_result) <- 'H3K9me3'
row.names(sig_H3K9me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.4715 odds ratio 1.071266 
#X-squared = 0.51464, df = 1, p-value = 0.4731






dmr_sig_H3K27ac_2 <- fread("dmr_sig_H3K27ac_2.bedGraph")
table(dmr_sig_H3K27ac_2$V20)
#Down-regulated   Up-regulated 
#  793           1593 
#Down-regulated   Up-regulated 
#    164            229 
sig_H3K27ac_2_fisher_data <- matrix(c(229, 1364,164, 629), nrow = 2)
colnames(sig_H3K27ac_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_2_fisher_data) <- c("H3K27ac_2", "other")
fisher.test(sig_H3K27ac_2_fisher_data) 
chisq.test(sig_H3K27ac_2_fisher_data)
sig_H3K27ac_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_2_fisher_data)))
sig_H3K27ac_2_result <- matrix(c(sig_H3K27ac_2_chisq_result[1:3,],sig_H3K27ac_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_2_result) <- 'H3K27ac_2'
row.names(sig_H3K27ac_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.0001349 odds ratio  0.6440207 
#X-squared = 14.846, df = 1, p-value = 0.0001167


dmr_sig_H3K27ac <- fread("dmr_sig_H3K27ac.bedGraph")
table(dmr_sig_H3K27ac$V14)
dmr_sig_H3K27ac_sub <- subset(dmr_sig_H3K27ac,subset=V4>0)
dmr_sig_H3K27ac_sub$region <- paste(dmr_sig_H3K27ac_sub$V5,dmr_sig_H3K27ac_sub$V6,sep="_")
dmr_sig_H3K27ac_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K27ac_sub$region))
table(dmr_sig_H3K27ac_sub$state)
#Down-regulated   Up-regulated 
#  793           1593 
#Down-regulated   Up-regulated 
#    641            1256
sig_H3K27ac_fisher_data <- matrix(c(1256, 337,641, 152), nrow = 2)
colnames(sig_H3K27ac_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_fisher_data) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac_fisher_data) 
chisq.test(sig_H3K27ac_fisher_data)
sig_H3K27ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_result <- matrix(c(sig_H3K27ac_chisq_result[1:3,],sig_H3K27ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_result) <- 'H3K27ac'
row.names(sig_H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 7.895e-07 odds ratio  0.6152394 
#X-squared = 24.858, df = 1, p-value = 6.171e-07



dmr_sig_H3K4me3 <- fread("dmr_sig_H3K4me3.bedGraph")
table(dmr_sig_H3K4me3$V14)
dmr_sig_H3K4me3_sub <- subset(dmr_sig_H3K4me3,subset=V4>3)
dmr_sig_H3K4me3_sub$region <- paste(dmr_sig_H3K4me3_sub$V5,dmr_sig_H3K4me3_sub$V6,sep="_")
dmr_sig_H3K4me3_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K4me3_sub$region))
table(dmr_sig_H3K4me3_sub$state)
#Down-regulated   Up-regulated 
#  793           1593 
#Down-regulated   Up-regulated 
#  337            557
sig_H3K4me3_fisher_data <- matrix(c(557, 1036,337 , 456), nrow = 2)
colnames(sig_H3K4me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_fisher_data) <- c("H3K4me3", "other")
fisher.test(sig_H3K4me3_fisher_data) 
chisq.test(sig_H3K4me3_fisher_data)
sig_H3K4me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_fisher_data)))
sig_H3K4me3_result <- matrix(c(sig_H3K4me3_chisq_result[1:3,],sig_H3K4me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_result) <- 'H3K4me3'
row.names(sig_H3K4me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.0003879 odds ratio  0.7276143 
#X-squared = 12.498, df = 1, p-value = 0.0004074


dmr_sig_H3K4me3_2 <- fread("dmr_sig_H3K4me3_2.bedGraph")
table(dmr_sig_H3K4me3_2$V14)
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig_H3K4me3_2,subset=V4>3)
dmr_sig_H3K4me3_2_sub$region <- paste(dmr_sig_H3K4me3_2_sub$V5,dmr_sig_H3K4me3_2_sub$V6,sep="_")
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K4me3_2_sub$region))
table(dmr_sig_H3K4me3_2_sub$state)
#Down-regulated   Up-regulated 
#  793           1593 
#Down-regulated   Up-regulated 
#   190            257 
sig_H3K4me3_2_fisher_data <- matrix(c(257, 1336,190, 603), nrow = 2)
colnames(sig_H3K4me3_2_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me3_2_fisher_data) <- c("H3K4me3_2", "other")
fisher.test(sig_H3K4me3_2_fisher_data) 
chisq.test(sig_H3K4me3_2_fisher_data)
sig_H3K4me3_2_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_2_fisher_data)))
sig_H3K4me3_2_result <- matrix(c(sig_H3K4me3_2_chisq_result[1:3,],sig_H3K4me3_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me3_2_result) <- 'H3K4me3_2'
row.names(sig_H3K4me3_2_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 6.078e-06 odds ratio  0.6106205 
#X-squared = 20.791, df = 1, p-value = 5.123e-06





###promoter
dmr_sig_promoter <- fread("dmr_sig_promoter.bedGraph")
dmr_sig_promoter$region <- paste(dmr_sig_promoter$V7,dmr_sig_promoter$V8,sep="_")
dmr_sig_promoter_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_promoter$region))
table(dmr_sig_promoter_sub$state)
#Down-regulated   Up-regulated 
#150            251 
#total
#Down-regulated   Up-regulated 
#  793           1593 

dmr_sig_promoter_fisher_data <- matrix(c(251, 1342, 150,643), nrow = 2)
colnames(dmr_sig_promoter_fisher_data) <- c("up", "down") 
rownames(dmr_sig_promoter_fisher_data) <- c("promoter", "other")
fisher.test(dmr_sig_promoter_fisher_data)
chisq.test(dmr_sig_promoter_fisher_data)
sig_promoter_chisq_result <- as.data.frame(unlist(chisq.test(dmr_sig_promoter_fisher_data)))
sig_promoter_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_promoter_fisher_data)))
sig_promoter_result <- matrix(c(sig_promoter_chisq_result[1:3,],sig_promoter_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_promoter_result) <- 'promoter'
row.names(sig_promoter_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.05516 odds ratio 0.8018472
#X-squared = 3.5564, df = 1, p-value = 0.05932




dmr_sig_H3K36me3 <- fread("dmr_sig_H3K36me3.bedGraph")
table(dmr_sig_H3K36me3$V20)
#Down-regulated   Up-regulated 
#  793           1593 
#Down-regulated   Up-regulated 
#    38            113 
sig_H3K36me3_fisher_data <- matrix(c(113, 1480,38, 655), nrow = 2)
colnames(sig_H3K36me3_fisher_data) <- c("up", "down") 
rownames(sig_H3K36me3_fisher_data) <- c("H3K36me3", "other")
fisher.test(sig_H3K36me3_fisher_data) 
chisq.test(sig_H3K36me3_fisher_data)
sig_H3K36me3_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K36me3_fisher_data)))
sig_H3K36me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K36me3_fisher_data)))
sig_H3K36me3_result <- matrix(c(sig_H3K36me3_chisq_result[1:3,],sig_H3K36me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K36me3_result) <- 'H3K36me3'
row.names(sig_H3K36me3_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.1696 odds ratio    1.315941 
#X-squared = 14.846, df = 1, p-value = 0.1825



dmr_sig_test_result <- as.data.frame(t(
  cbind(sig_promoter_result,
        sig_H3K27ac_2_result,sig_H3K27ac_result,
        sig_H3K9me3_result,
        sig_H3K4me3_result,sig_H3K4me3_2_result,sig_H3K36me3_result,
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

ggplot(dmr_sig_test_result[c('PMD','HMD'),],aes(x=lnR,y=logfdr))+
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























###total
table(dmrs$state)
#Down-regulated   Up-regulated 
#  14217          27255

###pmd-hmd
table(dmrs_pmd_sub_p$V16)
#Down-regulated   Up-regulated 
#  5239            9094 

table(dmrs_pmd_sub_h$V16)
#Down-regulated   Up-regulated 
# 3762            8580


sig_pmd_p_fisher_data <- matrix(c(454, 217, 1139, 576), nrow = 2)
colnames(sig_pmd_p_fisher_data) <- c("up", "down") 
rownames(sig_pmd_p_fisher_data) <- c("pmd", "other")
fisher.test(sig_pmd_p_fisher_data) 
chisq.test(sig_pmd_p_fisher_data)
sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_p_fisher_data) ))
#p-value = 0.595 odds ratio 1.058011 
#X-squared = 0.28374, df = 1, p-value = 0.5943


sig_pmd_h_fisher_data <- matrix(c(407, 170, 1186, 663), nrow = 2)
colnames(sig_pmd_h_fisher_data) <- c("up", "down") 
rownames(sig_pmd_h_fisher_data) <- c("hmd", "other")
fisher.test(sig_pmd_h_fisher_data) 
chisq.test(sig_pmd_h_fisher_data)
sig_pmd_h_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_h_fisher_data) ))
#p-value = 0.004905  odds ratio 1.338205 
#X-squared = 7.6943, df = 1, p-value = 0.005539



sig_pmd_pvh_fisher_data <- matrix(c(454, 217, 407, 170), nrow = 2)
colnames(sig_pmd_pvh_fisher_data) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")

mcnemar.test(sig_pmd_pvh_fisher_data, correct=F) 

sig_pmd_pvh_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_pvh_fisher_data) ))
#p-value = 0.2968  odds ratio 0.8739731 
#McNemar's chi-squared = 57.853, df = 1, p-value = 2.825e-14

sig_pmd_pvh_fisher_data2 <- matrix(c(9094,8580,5239, 3762), nrow = 2)
colnames(sig_pmd_pvh_fisher_data2) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data2) <- c("pmd", "hmd")

mcnemar.test(sig_pmd_pvh_fisher_data2, correct=F) 
#McNemar's chi-squared = 807.75, df = 1, p-value < 2.2e-16 odds_ratio=0.7610933




















###total
table(dmr_sig$state)
#Down-regulated   Up-regulated 
#  700           1338  

###pmd-hmd
table(dmr_sig_pmd_sub_p$V16)
#Down-regulated   Up-regulated 
# 194            430 

table(dmr_sig_pmd_sub_h$V16)
#Down-regulated   Up-regulated 
#146            318 


test_fisher_data <- matrix(c(194,506,430, 908), nrow = 2)
colnames(test_fisher_data) <- c("down", "up") 
rownames(test_fisher_data) <- c("pmd", "other")
fisher.test(test_fisher_data) 
chisq.test(test_fisher_data)

table(dmr_pmd_sub_p$V16)
#Down-regulated   Up-regulated 
#  217            454 

table(dmr_pmd_sub_h$V16)
#Down-regulated   Up-regulated 
# 170            407


test_fisher_data2 <- matrix(c(454, 217, 407, 170), nrow = 2)
colnames(sig_pmd_p_fisher_data) <- c("up", "down") 
rownames(sig_pmd_p_fisher_data) <- c("pmd", "other")
fisher.test(sig_pmd_p_fisher_data) 
chisq.test(sig_pmd_p_fisher_data)
sig_pmd_p_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_p_fisher_data)))
sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_p_fisher_data)))
sig_pmd_p_result <- matrix(c(sig_pmd_p_chisq_result[1:3,],sig_pmd_p_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_p_result) <- 'PMD'
row.names(sig_pmd_p_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.595 odds ratio 1.058011 
#X-squared = 0.28374, df = 1, p-value = 0.5943





test_h3k27ac <- fread("H3K27ac_PMD.bedGraph")
test_h3k9me3 <- fread("H3K9me3_PMD.bedGraph")
test_h3k4me3 <- fread("H3K4me3_PMD.bedGraph")
test_h3k36me3 <- fread("H3K36me3_PMD.bedGraph")


PMD_2 <- dplyr::select(PMD,1,2,3,6)
write.table(PMD_2, "hg38_PMD.bed",sep = "\t",quote=F,row.names = F,col.names = F)
dmr_sig_up <- subset(dmr_sig,subset=state == "Up-regulated")
dmr_sig_up <-dplyr::select(dmr_sig_up,1,2,3,10)
dmr_sig_down <- subset(dmr_sig,subset=state == "Down-regulated")
dmr_sig_down <-dplyr::select(dmr_sig_down,1,2,3,10)
write.table(dmr_sig_up, "dmr_sig_up.bed",sep = "\t",quote=F,row.names = F,col.names = F)
write.table(dmr_sig_down, "dmr_sig_down.bed",sep = "\t",quote=F,row.names = F,col.names = F)


h3k27ac_down <- subset(dmr_sig_H3K27ac,subset=V14 == "Down-regulated")
h3k27ac_up <- subset(dmr_sig_H3K27ac,subset=V14 == "Up-regulated")

table(h3k27ac_up$V14)
table(h3k27ac_down$V14)

sum(h3k27ac_up$V4)
sum(h3k27ac_down$V4)
