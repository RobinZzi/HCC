library(DSS)
library(ggrepel)
library(dplyr)
library(stringr)
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_cov")
rm(list = ls())
load("dss.Rdata")
save.image("dss.Rdata")

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_cov")
GA45_c1 <- fread("GA45_V1.bismark.cov.gz")
GA45_c2 <- fread("GA45_V2.bismark.cov.gz")
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

dmls.sm <- callDML(dmlTest, p.threshold=0.05)
dmls.sm <- mutate(dmls.sm,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                            diff>0 ~ "Down-regulated", # 下调
                                            TRUE ~ "Unchanged"))


dmrs.sm <- callDMR(dmlTest.sm, p.threshold=0.05)
dmrs.sm <- mutate(dmrs.sm,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                            diff.Methy>0 ~ "Down-regulated", # 下调
                                            TRUE ~ "Unchanged"))
table(dmrs.sm$state)

dmrs.sm_up <- subset(dmrs.sm,subset=diff.Methy<0)
dmrs.sm_down <- subset(dmrs.sm,subset=diff.Methy>0)


dmlTest.sm_sig <- subset(dmlTest.sm,subset=fdr<0.05)
dmr_sig <- callDMR(dmlTest.sm_sig)
dmr_sig <- mutate(dmlTest.sm_sig_dmr,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                                                  diff.Methy>0 ~ "Down-regulated", # 下调
                                                                  TRUE ~ "Unchanged"))

dmr_sig$region <- paste(dmr_sig$chr,dmr_sig$start,sep = "_")

dmls_sig <- callDML(dmlTest.sm_sig)
dmls_sig <- mutate(dmls_sig,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                              diff>0 ~ "Down-regulated", # 下调
                                              TRUE ~ "Unchanged"))
table(dmls_sig$state)
showOneDMR(dmrs[1,], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['2408',], BSobj)


test <- fread("GA45_dmrs_sig.bedGraph")

write.table(dmr_sig, "GA45_dmrs_sig.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

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



dmrs.sm$region <- paste(dmrs.sm$chr,dmrs.sm$start,sep = "_")

dmr_pmd <- fread("dmr_pmd.bedGraph")


dmr_H3K9me3 <- fread("dmr_H3K9me3.bedGraph")


dmr_H3K27ac <- fread("dmr_H3K27ac.bedGraph")


dmr_H3K36me3 <- fread("dmr_H3K36me3.bedGraph")


dmr_H3K4me3 <- fread("dmr_H3K4me3.bedGraph")


dmr_promoter <- fread("dmr_promoter.bedGraph")



dmr_H3K4me3_2 <- fread("dmr_H3K4me3_2.bedGraph")
dmr_H3K4me3_2_sub <- subset(dmr_H3K4me3_2,subset=V4>3)
dmr_H3K4me3_2_sub$region <- paste(dmr_H3K4me3_2_sub$V5,dmr_H3K4me3_2_sub$V6,sep="_")
dmr_H3K4me3_2_sub <- subset(dmrs.sm,subset=region %in% unique(dmr_H3K4me3_2_sub$region))


dmr_H3K27ac_2 <- fread("dmr_H3K27ac_2.bedGraph")
dmr_H3K27ac_2 <- fread("dmr_H3K27ac_2.bedGraph")
dmr_H3K27ac_2_sub <- subset(dmr_H3K27ac_2,subset=V4>3)
dmr_H3K27ac_2_sub$region <- paste(dmr_H3K27ac_2_sub$V5,dmr_H3K27ac_2_sub$V6,sep="_")
dmr_H3K27ac_2_sub <- subset(dmrs.sm,subset=region %in% unique(dmr_H3K27ac_2_sub$region))




dmr_pmd_sub <- subset(dmr_pmd,subset=V6!="Neither")
dmr_pmd_sub_up <- subset(dmr_pmd_sub,subset=V16=="Up-regulated")
dmr_pmd_sub_down <- subset(dmr_pmd_sub,subset=V16=="Down-regulated")


dmr_pmd_sub_p <- subset(dmr_pmd_sub,subset=V6=="commonPMD")

dmr_pmd_sub_h <- subset(dmr_pmd_sub,subset=V6=="commonHMD")



###total
table(dmrs.sm$state)
#Down-regulated   Up-regulated 
#46191          77597 

###pmd-hmd
table(dmr_pmd_sub_p$V16)
#Down-regulated   Up-regulated 
# 17229          29283

table(dmr_pmd_sub_h$V16)
#Down-regulated   Up-regulated 
#13043          22395 


sig_pmd_p_fisherdata <- matrix(c(17229, 29283, 28962, 48314), nrow = 2)
colnames(sig_pmd_p_fisherdata) <- c("up", "down") 
rownames(sig_pmd_p_fisherdata) <- c("pmd", "other")
fisher.test(sig_pmd_p_fisherdata) 
sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_p_fisherdata) ))
#p-value = 0.1248 odds ratio 0.9814943 

sig_pmd_h_fisherdata <- matrix(c(13043, 22395, 33148, 55202), nrow = 2)
colnames(sig_pmd_h_fisherdata) <- c("up", "down") 
rownames(sig_pmd_h_fisherdata) <- c("hmd", "other")
fisher.test(sig_pmd_h_fisherdata)  
sig_pmd_h_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_h_fisherdata) ))
#p-value = 0.01894  odds ratio 0.9698883 


###H3K9me3
table(dmr_H3K9me3$V20)
#Down-regulated   Up-regulated 
#33             67
#total
#Down-regulated   Up-regulated 
#46191          77597 

sig_H3K9me3_fisherdata <- matrix(c(33, 67, 46185, 77530), nrow = 2)
colnames(sig_H3K9me3_fisherdata) <- c("up", "down") 
rownames(sig_H3K9me3_fisherdata) <- c("H3K9me3", "other")
fisher.test(sig_H3K9me3_fisherdata) 
sig_H3K9me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K9me3_fisherdata) ))
#p-value = 0.4088 odds ratio 0.8268457




###H3K27ac
table(dmr_H3K27ac$V20)
#Down-regulated   Up-regulated 
#2907           4845
#total
#Down-regulated   Up-regulated 
#46191          77597  

sig_H3K27ac_fisherdata <- matrix(c(2907, 4845, 43284, 72752), nrow = 2)
colnames(sig_H3K27ac_fisherdata) <- c("up", "down") 
rownames(sig_H3K27ac_fisherdata) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac_fisherdata) 
sig_H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_fisherdata) ))
#p-value = 0.7342 odds ratio 1.008485



###H3K27ac_2
table(dmr_H3K27ac_2_sub$state)
#Down-regulated   Up-regulated 
#4734           7758
#total
#Down-regulated   Up-regulated 
#46191          77597  

sig_H3K27ac_2_fisherdata <- matrix(c(4734, 7758, 41457, 69839), nrow = 2)
colnames(sig_H3K27ac_2_fisherdata) <- c("up", "down") 
rownames(sig_H3K27ac_2_fisherdata) <- c("H3K27ac_2", "other")
fisher.test(sig_H3K27ac_2_fisherdata) 
sig_H3K27ac_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_2_fisherdata) ))
#p-value = 0.1572 odds ratio 1.027983





###H3K36me3
table(dmr_H3K36me3$V20)
#Down-regulated   Up-regulated 
#2824           4882
#total
#Down-regulated   Up-regulated 
#46191          77597  

sig_H3K36me3_fisherdata <- matrix(c(2824, 4882, 43367, 72715), nrow = 2)
colnames(sig_H3K36me3_fisherdata) <- c("up", "down") 
rownames(sig_H3K36me3_fisherdata) <- c("H3K36me3", "other")
fisher.test(sig_H3K36me3_fisherdata) 
sig_H3K36me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K36me3_fisherdata) ))
#p-value = 0.2148 odds ratio 0.9699085




###H3K4me3
table(dmr_H3K4me3$V20)
#Down-regulated   Up-regulated 
#3627           5708
#total
#Down-regulated   Up-regulated 
#46191          77597  

sig_H3K4me3_fisherdata <- matrix(c(3627, 5708, 42564, 71889), nrow = 2)
colnames(sig_H3K4me3_fisherdata) <- c("up", "down") 
rownames(sig_H3K4me3_fisherdata) <- c("H3K4me3", "other")
fisher.test(sig_H3K4me3_fisherdata) 
sig_H3K4me3_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_fisherdata) ))
#p-value = 0.001458 odds ratio 1.073206

###H3K4me3_2
table(dmr_H3K4me3_2_sub$state)
#Down-regulated   Up-regulated 
#3218           5165
#total
#Down-regulated   Up-regulated 
#46191          77597  

sig_H3K4me3_2_fisherdata <- matrix(c(3218, 5165, 42973, 72432), nrow = 2)
colnames(sig_H3K4me3_2_fisherdata) <- c("up", "down") 
rownames(sig_H3K4me3_2_fisherdata) <- c("H3K4me3_2", "other")
fisher.test(sig_H3K4me3_2_fisherdata) 
sig_H3K4me3_2_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me3_2_fisherdata) ))
#p-value = 0.03632 odds ratio 1.050135



###promoter
table(dmr_promoter$V16)
#Down-regulated   Up-regulated 
#5621           8831
#total
#Down-regulated   Up-regulated 
#46191          77597  

sig_promoter_fisherdata <- matrix(c(5621, 8831, 40570, 68766), nrow = 2)
colnames(sig_promoter_fisherdata) <- c("up", "down") 
rownames(sig_promoter_fisherdata) <- c("promoter", "other")
sig_promoter_fisher_result <- as.data.frame(unlist(fisher.test(sig_promoter_fisherdata))) 
#p-value = 3.131e-05 odds ratio 1.078877








fisher_test_result <- as.data.frame(t(
  cbind(sig_pmd_p_fisher_result,sig_pmd_h_fisher_result,sig_H3K9me3_fisher_result,sig_H3K27ac_fisher_result,sig_H3K27ac_2_fisher_result,
        sig_H3K4me3_2_fisher_result,sig_H3K4me3_fisherdata_result,sig_H3K36me3_fisher_result,sig_promoter_fisher_result)))
fisher_test_result$histone_mark <- c('pmd','hmd','H3K9me3','H3K27ac','H3K27ac_2','H3K4me3_2','H3K4me3','H3K36me3','promoter')
row.names(fisher_test_result) <- fisher_test_result$histone_mark
fisher_test_result <- dplyr::select(fisher_test_result,1,4,9)

colnames(fisher_test_result) <- c("P_value","odds_ratio","Histone_Mark")
fisher_test_result$P_value <- as.numeric(fisher_test_result$P_value)
fisher_test_result$odds_ratio <- as.numeric(fisher_test_result$odds_ratio)

fisher_test_result[,"logp"] <- -log10(fisher_test_result$P_value)
fisher_test_result[,"lnR"] <- -log(fisher_test_result$odds_ratio)

fisher_test_result <- dplyr::mutate(fisher_test_result,case = case_when(lnR > 0 ~ 'Enriched in Hyper',
                                                                        lnR < 0 ~ 'Enriched in Hypo'))



ggplot(fisher_test_result,aes(x=lnR,y=logp))+
  geom_point(size = 1,aes(color = case))+# 根据expression水平进行着色
  xlab(expression("ln odds ratio")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" P")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-1, 1)) +
  theme_bw() +
  geom_text_repel(aes(label=Histone_Mark, color = case), size =5,hjust=1,vjust=1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Fisher test result of Histone Mark ")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")






















##
dmr_sig_pmd <- fread("dmr_sig_pmd.bedGraph")


dmr_sig_H3K9me3 <- fread("dmr_sig_H3K9me3.bedGraph")

dmr_sig_H3K9me3_2 <- fread("dmr_sig_H3K9me3_2.bedGraph")
dmr_sig_H3K9me3_2_sub <- subset(dmr_sig_H3K9me3_2,subset=V4>3)
dmr_sig_H3K9me3_2_sub$region <- paste(dmr_sig_H3K9me3_2_sub$V5,dmr_sig_H3K9me3_2_sub$V6,sep="_")
dmr_sig_H3K9me3_2_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K9me3_2_sub$region))

dmr_sig_H3K9me3_3 <- fread("dmr_sig_H3K9me3_3.bedGraph")

dmr_sig_H3K27ac <- fread("dmr_sig_H3K27ac.bedGraph")


dmr_sig_H3K36me3 <- fread("dmr_sig_H3K36me3.bedGraph")


dmr_sig_H3K4me3 <- fread("dmr_sig_H3K4me3.bedGraph")


dmr_sig_promoter <- fread("dmr_sig_promoter.bedGraph")
dmr_sig_promoter$region <- paste(dmr_sig_promoter$V7,dmr_sig_promoter$V8,sep="_")
dmr_sig_promoter_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_promoter$region))


dmr_sig_H3K4me3_2 <- fread("dmr_sig_H3K4me3_2.bedGraph")
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig_H3K4me3_2,subset=V4>3)
dmr_sig_H3K4me3_2_sub$region <- paste(dmr_sig_H3K4me3_2_sub$V5,dmr_sig_H3K4me3_2_sub$V6,sep="_")
dmr_sig_H3K4me3_2_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K4me3_2_sub$region))


dmr_sig_H3K4me3_3 <- fread("dmr_sig_H3K4me3_3.bedGraph")
dmr_sig_H3K4me3_3_sub <- subset(dmr_sig_H3K4me3_3,subset=V4>3)
dmr_sig_H3K4me3_3_sub$region <- paste(dmr_sig_H3K4me3_3_sub$V5,dmr_sig_H3K4me3_3_sub$V6,sep="_")
dmr_sig_H3K4me3_3_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K4me3_3_sub$region))




dmr_sig_H3K27ac_2 <- fread("dmr_sig_H3K27ac_2.bedGraph")
dmr_sig_H3K27ac_2 <- fread("dmr_sig_H3K27ac_2.bedGraph")
dmr_sig_H3K27ac_2_sub <- subset(dmr_sig_H3K27ac_2,subset=V4>3)
dmr_sig_H3K27ac_2_sub$region <- paste(dmr_sig_H3K27ac_2_sub$V5,dmr_sig_H3K27ac_2_sub$V6,sep="_")
dmr_sig_H3K27ac_2_sub <- subset(dmr_sig,subset=region %in% unique(dmr_sig_H3K27ac_2_sub$region))




dmr_sig_pmd_sub <- subset(dmr_sig_pmd,subset=V6!="Neither")
dmr_sig_pmd_sub_up <- subset(dmr_sig_pmd_sub,subset=V16=="Up-regulated")
dmr_sig_pmd_sub_down <- subset(dmr_sig_pmd_sub,subset=V16=="Down-regulated")


dmr_sig_pmd_sub_p <- subset(dmr_sig_pmd_sub,subset=V6=="commonPMD")

dmr_sig_pmd_sub_h <- subset(dmr_sig_pmd_sub,subset=V6=="commonHMD")



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


dmr_sig_pmd_p_fisher_data <- matrix(c(194, 430, 506, 908), nrow = 2)
colnames(dmr_sig_pmd_p_fisher_data) <- c("down", "up") 
rownames(dmr_sig_pmd_p_fisher_data) <- c("pmd", "other")
fisher.test(dmr_sig_pmd_p_fisher_data) 
dmr_sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_pmd_p_fisherdata) ))
#p-value = 0.04292 odds ratio 0.8096806 

dmr_sig_pmd_h_fisher_data <- matrix(c(146, 318, 554, 1020), nrow = 2)
colnames(dmr_sig_pmd_h_fisher_data) <- c("down", "up") 
rownames(dmr_sig_pmd_h_fisher_data) <- c("hmd", "other")
fisher.test(dmr_sig_pmd_h_fisherdata)  
 dmr_sig_pmd_h_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_pmd_h_fisher_data) ))
#p-value = 0.1481  odds ratio 0.8453794 

dmr_sig_pmd_pvh_fisher_data <- matrix(c(194, 430, 146, 318), nrow = 2)
colnames(dmr_sig_pmd_pvh_fisher_data) <- c("down", "up") 
rownames(dmr_sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")
fisher.test(dmr_sig_pmd_pvh_fisher_data) 
dmr_sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_pmd_p_fisherdata) ))
#p-value = 0.04292 odds ratio 0.8096806 





###H3K9me3
table(dmr_sig_H3K9me3$V20)
#Down-regulated   Up-regulated 
#3             4
#total
#Down-regulated   Up-regulated 
#700           1338 

dmr_sig_H3K9me3_fisher_data <- matrix(c(15, 20, 685, 1318), nrow = 2)
colnames(dmr_sig_H3K9me3_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K9me3_fisher_data) <- c("H3K9me3", "other")
fisher.test(dmr_sig_H3K9me3_fisher_data) 
dmr_sig_H3K9me3_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K9me3_fisher_data) ))
#p-value = 0.2863 odds ratio 1.442847



###H3K9me3_2
table(dmr_sig_H3K9me3_2_sub$state)
#Down-regulated   Up-regulated 
#259            486
#total
#Down-regulated   Up-regulated 
#700           1338 

dmr_sig_H3K9me3_2_fisher_data <- matrix(c(259, 486, 441, 852), nrow = 2)
colnames(dmr_sig_H3K9me3_2_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K9me3_2_fisher_data) <- c("H3K9me3", "other")
fisher.test(dmr_sig_H3K9me3_2_fisher_data) 
dmr_sig_H3K9me3_2_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K9me3_2_fisher_data) ))
#p-value = 0.6967 odds ratio 1.443502



###H3K27ac
table(dmr_sig_H3K27ac$V20)
#Down-regulated   Up-regulated 
#119           251
#total
#Down-regulated   Up-regulated 
#700           1338  

dmr_sig_H3K27ac_fisher_data <- matrix(c(119, 251, 581, 1087), nrow = 2)
colnames(dmr_sig_H3K27ac_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K27ac_fisher_data) <- c("H3K27ac", "other")
fisher.test(dmr_sig_H3K27ac_fisher_data) 
dmr_sig_H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K27ac_fisher_data) ))
#p-value = 0.3338 odds ratio 0.887028



###H3K27ac_2
table(dmr_sig_H3K27ac_2_sub$state)
#Down-regulated   Up-regulated 
#184           360
#total
#Down-regulated   Up-regulated 
#700           1338  

dmr_sig_H3K27ac_2_fisher_data <- matrix(c(184, 360, 516, 978), nrow = 2)
colnames(dmr_sig_H3K27ac_2_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K27ac_2_fisher_data) <- c("H3K27ac_2", "other")
fisher.test(dmr_sig_H3K27ac_2_fisher_data) 
dmr_sig_H3K27ac_2_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K27ac_2_fisher_data) ))
#p-value = 0.7921 odds ratio 0.9687447





###H3K36me3
table(dmr_sig_H3K36me3$V20)
#Down-regulated   Up-regulated 
#33           86
#total
#Down-regulated   Up-regulated 
#700           1338  

dmr_sig_H3K36me3_fisher_data <- matrix(c(33, 86, 667, 1252), nrow = 2)
colnames(dmr_sig_H3K36me3_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K36me3_fisher_data) <- c("H3K36me3", "other")
fisher.test(dmr_sig_H3K36me3_fisher_data) 
dmr_sig_H3K36me3_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K36me3_fisher_data) ))
#p-value = 0.1355 odds ratio 0.7203782




###H3K4me3
table(dmr_sig_H3K4me3$V20)
#Down-regulated   Up-regulated 
#157           312
#total
#Down-regulated   Up-regulated 
#700           1338   

dmr_sig_H3K4me3_fisher_data <- matrix(c(157, 312, 543, 1026), nrow = 2)
colnames(dmr_sig_H3K4me3_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K4me3_fisher_data) <- c("H3K4me3", "other")
fisher.test(dmr_sig_H3K4me3_fisher_data) 
dmr_sig_H3K4me3_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K4me3_fisher_data) ))
#p-value = 0.658 odds ratio 0.9508225

###H3K4me3_2
table(dmr_sig_H3K4me3_2_sub$state)
#Down-regulated   Up-regulated 
#133           259
#total
#Down-regulated   Up-regulated 
#700          1338  

dmr_sig_H3K4me3_2_fisher_data <- matrix(c(133, 259, 567, 1079), nrow = 2)
colnames(dmr_sig_H3K4me3_2_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K4me3_2_fisher_data) <- c("H3K4me3_2", "other")
fisher.test(dmr_sig_H3K4me3_2_fisher_data) 
dmr_sig_H3K4me3_2_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K4me3_2_fisher_data) ))
#p-value = 0.8593 odds ratio 0.9772248


###H3K4me3_3
table(dmr_sig_H3K4me3_3_sub$state)
#Down-regulated   Up-regulated 
#279           523
#total
#Down-regulated   Up-regulated 
#700          1338  

dmr_sig_H3K4me3_3_fisher_data <- matrix(c(279, 523, 421, 815), nrow = 2)
colnames(dmr_sig_H3K4me3_3_fisher_data) <- c("down", "up") 
rownames(dmr_sig_H3K4me3_3_fisher_data) <- c("H3K4me3_2", "other")
fisher.test(dmr_sig_H3K4me3_3_fisher_data) 
dmr_sig_H3K4me3_3_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_H3K4me3_3_fisher_data) ))
#p-value = 0.7385 odds ratio 1.032713


###promoter
table(dmr_sig_promoter_sub$state)
#Down-regulated   Up-regulated 
#119           253
#total
#Down-regulated   Up-regulated 
#700          1338 

dmr_sig_promoter_fisher_data <- matrix(c(119, 253, 581, 1085), nrow = 2)
colnames(dmr_sig_promoter_fisher_data) <- c("up", "down") 
rownames(dmr_sig_promoter_fisher_data) <- c("promoter", "other")
fisher.test(dmr_sig_promoter_fisher_data)
dmr_sig_promoter_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_promoter_fisher_data))) 
#p-value = 0.3049 odds ratio 0.8784575







dmr_sig_fisher_test_result <- as.data.frame(t(
  cbind(dmr_sig_pmd_p_fisher_result,dmr_sig_pmd_h_fisher_result,dmr_sig_H3K9me3_fisher_result,
        dmr_sig_H3K27ac_fisher_result,dmr_sig_H3K27ac_2_fisher_result,dmr_sig_H3K4me3_2_fisher_result,
        dmr_sig_H3K4me3_fisher_result,dmr_sig_H3K36me3_fisher_result,dmr_sig_promoter_fisher_result,
        dmr_sig_H3K4me3_3_fisher_result,dmr_sig_H3K9me3_2_fisher_result)))
dmr_sig_fisher_test_result$histone_mark <- c('pmd','hmd','H3K9me3',
                                             'H3K27ac','H3K27ac_2','H3K4me3_2',
                                             'H3K4me3','H3K36me3','promoter',
                                             'H3K4me3_3','H3K9me3_2')
row.names(dmr_sig_fisher_test_result) <- dmr_sig_fisher_test_result$histone_mark
dmr_sig_fisher_test_result <- dplyr::select(dmr_sig_fisher_test_result,1,4,9)

colnames(dmr_sig_fisher_test_result) <- c("P_value","odds_ratio","Histone_Mark")
dmr_sig_fisher_test_result$P_value <- as.numeric(dmr_sig_fisher_test_result$P_value)
dmr_sig_fisher_test_result$odds_ratio <- as.numeric(dmr_sig_fisher_test_result$odds_ratio)

dmr_sig_fisher_test_result[,"logp"] <- -log10(dmr_sig_fisher_test_result$P_value)
dmr_sig_fisher_test_result[,"lnR"] <- -log(dmr_sig_fisher_test_result$odds_ratio)

dmr_sig_fisher_test_result <- dplyr::mutate(dmr_sig_fisher_test_result,case = case_when(lnR > 0 ~ 'Enriched in Hyper',
                                                                                        lnR < 0 ~ 'Enriched in Hypo'))

dmr_sig_fisher_test_result$fdr <- p.adjust(dmr_sig_fisher_test_result$P_value, method  = "BH")


ggplot(dmr_sig_fisher_test_result,aes(x=lnR,y=logp))+
  geom_point(size = 1,aes(color = case))+# 根据expression水平进行着色
  xlab(expression("ln odds ratio")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" P")) + # 修饰y轴题目
  theme_bw() +
  geom_text_repel(aes(label=Histone_Mark, color = case), size =5,hjust=1,vjust=1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Fisher test result of Histone Mark ")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")









h3k9me3_bed <- fread("hg38_293T_H3K9me3.bedGraph")
h3k9me3_bed <- subset(h3k9me3_bed,subset = V1 != 'NC_012920.1')
h3k9me3_bed$chr <- str_extract(h3k9me3_bed$V1, "\\d{2}(?=\\.)")
h3k9me3_bed$chr <- as.numeric(h3k9me3_bed$chr)
h3k9me3_bed$chr <- paste("chr",h3k9me3_bed$chr,sep = "")
h3k9me3_bed$chr <- gsub("chr23", "chrX", h3k9me3_bed$chr)
h3k9me3_bed$chr <- gsub("chr24", "chrY", h3k9me3_bed$chr)
h3k9me3_bed <- dplyr::select(h3k9me3_bed,5,2,3,4)
colnames(h3k9me3_bed) <- c('chr', 'start','end','level')
write.table(h3k9me3_bed, "293T_H3K9me3.bedgraph",sep = "\t",quote=F,row.names = F,col.names = F)




hek4me3_bed <- fread("hg38_293T_H3K4me3_2.bedGraph")
hek4me3_bed <- subset(hek4me3_bed,subset = V1 != 'NC_012920.1')
hek4me3_bed$chr <- str_extract(hek4me3_bed$V1, "\\d{2}(?=\\.)")
hek4me3_bed$chr <- as.numeric(hek4me3_bed$chr)
hek4me3_bed$chr <- paste("chr",hek4me3_bed$chr,sep = "")
hek4me3_bed$chr <- gsub("chr23", "chrX", hek4me3_bed$chr)
hek4me3_bed$chr <- gsub("chr24", "chrY", hek4me3_bed$chr)
hek4me3_bed <- dplyr::select(hek4me3_bed,5,2,3,4)
colnames(hek4me3_bed) <- c('chr', 'start','end','level')
write.table(hek4me3_bed, "293T_H3K4me3.bedgraph",sep = "\t",quote=F,row.names = F,col.names = F)







PMD <- fread("PMD_coordinates_hg38.bed")
