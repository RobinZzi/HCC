library(DSS)
library(ggrepel)
library(dplyr)

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_cov")
rm(list = ls())
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
dmlTest.sm_sig_dmr <- callDMR(dmlTest.sm_sig)
dmlTest.sm_sig_dmr <- mutate(dmlTest.sm_sig_dmr,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                                                  diff.Methy>0 ~ "Down-regulated", # 下调
                                                                  TRUE ~ "Unchanged"))



dmls_sig <- callDML(dmlTest.sm_sig)
dmls_sig <- mutate(dmls_sig,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                              diff>0 ~ "Down-regulated", # 下调
                                              TRUE ~ "Unchanged"))
table(dmls_sig$state)
showOneDMR(dmrs[1,], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['1179',], BSobj)






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


dmls_sig_prop <- as.data.frame(t(table(dmlTest.sm_sig_dmr$state)))
ggplot(dmls_sig_prop, aes(x = 1, y = Freq, fill = Var2))+
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




dmr_pmd <- fread("dmr_pmd.bedGraph")


dmr_H3K9me3 <- fread("dmr_H3K9me3.bedGraph")


dmr_H3K27ac <- fread("dmr_H3K27ac.bedGraph")


dmr_H3K36me3 <- fread("dmr_H3K36me3.bedGraph")


dmr_H3K4me3 <- fread("dmr_H3K4me3.bedGraph")


dmr_promoter <- fread("dmr_promoter.bedGraph")




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
#p-value = 0.1248 odds ratio 0.9814943 

sig_pmd_h_fisherdata <- matrix(c(13043, 22395, 33148, 55202), nrow = 2)
colnames(sig_pmd_h_fisherdata) <- c("up", "down") 
rownames(sig_pmd_h_fisherdata) <- c("hmd", "other")
fisher.test(sig_pmd_h_fisherdata)  
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
#p-value = 0.7342 odds ratio 1.008485




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
#p-value = 0.001458 odds ratio 1.073206




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
fisher.test(sig_promoter_fisherdata) 
#p-value = 3.131e-05 odds ratio 1.078877



sh_v_down <- subset(sh_v_diff,subset=sh_v_diff$logFC < 0)
sh_v_up <- subset(sh_v_diff,subset=sh_v_diff$logFC > 0)

pro_meth_up <- subset(dmr_sig_promoter,subset=V14 < 0)
pro_meth_down <- subset(dmr_sig_promoter,subset=V14 > 0)



common_up <- intersect(sh_v_up$GeneID,unique(pro_meth_down$V5))
common_down <- intersect(sh_v_down$GeneID,unique(pro_meth_up$V5))






fisher_test_result <- as.data.frame(
  cbind(c("PMD","HMD","H3K9me3","H3K27me3","H3K36me3","H3K27ac-1","H3K27ac-2","H3K4me1","H3K9ac","H3K36ac","H4K5ac","promoter"),
        c(0.0504,0.003241,0.9047,0.2304,1,1.545e-12,4.472e-10,0.1548,1.776e-09,7.969e-09,6.066e-09,0.0018),
        c(1.173642,0.7886059,0.9574356,1.267434,0.988317,0.5346163,0.6033645,1.133057,0.6071659,0.626763,0.6204754,0.7768482)
  )
)

colnames(fisher_test_result) <- c("Histone_Mark","P_value","odds_ratio")
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
  geom_text_repel(aes(label=Histone_Mark, color = case), size =3,hjust=1,vjust=-1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Fisher test result of Hitone ")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") 



