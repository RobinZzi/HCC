library(DSS)
library(ggrepel)
library(dplyr)

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240611/methy_mtx")
rm(list = ls())

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240611/meth_ex")
SNHG6_c1 <- fread("SNHG6_V1_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
SNHG6_c2 <- fread("SNHG6_V2_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
SNHG6_e1 <- fread("SNHG6_SH2-1_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240614/meth_ex")
SNHG6_e2 <- fread("SNHG6_E2_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")

save.image("dss.Rdata")

SNHG6_c1$total <- SNHG6_c1$V5+SNHG6_c1$V6
SNHG6_c1_con <- dplyr::select(SNHG6_c1,1,2,7,5)
colnames(SNHG6_c1_con) <- c("chr","pos","N","X")

SNHG6_c2$total <- SNHG6_c2$V5+SNHG6_c2$V6
SNHG6_c2_con <- dplyr::select(SNHG6_c2,1,2,7,5)
colnames(SNHG6_c2_con) <- c("chr","pos","N","X")

SNHG6_e1$total <- SNHG6_e1$V5+SNHG6_e1$V6
SNHG6_e1_con <- dplyr::select(SNHG6_e1,1,2,7,5)
colnames(SNHG6_e1_con) <- c("chr","pos","N","X")

SNHG6_e2$total <- SNHG6_e2$V5+SNHG6_e2$V6
SNHG6_e2_con <- dplyr::select(SNHG6_e2,1,2,7,5)
colnames(SNHG6_e2_con) <- c("chr","pos","N","X")

BSobj <- makeBSseqData( list(SNHG6_c1_con, SNHG6_c2_con, SNHG6_e1_con, SNHG6_e2_con),
                        c("C1","C2", "E1", "E2") )
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("E1", "E2"))
dmlTest.sm <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("E1", "E2"), smoothing=TRUE)

dmls.sm <- callDML(dmlTest, p.threshold=0.05)
dmls.sm <- mutate(dmls.sm,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                            diff>0 ~ "Down-regulated", # 下调
                                            TRUE ~ "Unchanged"))

dmls.sm2 <- callDML(dmlTest,delta=0.1, p.threshold=0.05)
dmls.sm2 <- mutate(dmls.sm2,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                            diff>0 ~ "Down-regulated", # 下调
                                            TRUE ~ "Unchanged"))


dmrs.sm <- callDMR(dmlTest.sm, p.threshold=0.01)
dmrs.sm <- mutate(dmrs.sm,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                            diff.Methy>0 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged"))
dmrs.sm_up <- subset(dmrs.sm,subset=diff.Methy<0)
dmrs.sm_down <- subset(dmrs.sm,subset=diff.Methy>0)

dmrs.sm2 <- callDMR(dmlTest.sm, delta=0.1, p.threshold=0.05)
dmrs.sm2 <- mutate(dmrs.sm2,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                            diff.Methy>0 ~ "Down-regulated", # 下调
                                            TRUE ~ "Unchanged"))

dmlTest.sm_sig <- subset(dmlTest.sm,subset=fdr<0.05)
dmlTest.sm_sig_dmr <- callDMR(dmlTest.sm_sig)
dmlTest.sm_sig_dmr <- mutate(dmlTest.sm_sig_dmr,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                              diff.Methy>0 ~ "Down-regulated", # 下调
                                              TRUE ~ "Unchanged"))

dmlTest.sm_sig2 <- subset(dmlTest.sm,subset=fdr<0.05)
dmlTest.sm_sig2_dmr <- callDMR(dmlTest.sm_sig2,delta=0.1)
dmlTest.sm_sig2_dmr <- mutate(dmlTest.sm_sig2_dmr,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                                                  diff.Methy>0 ~ "Down-regulated", # 下调
                                                                  TRUE ~ "Unchanged"))


dmls_sig <- callDML(dmlTest.sm_sig)
dmls_sig <- mutate(dmls_sig,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                            diff>0 ~ "Down-regulated", # 下调
                                            TRUE ~ "Unchanged"))

showOneDMR(dmrs[1,], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['1739',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['3166',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['3219',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['234',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['1739',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['2672',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['1417',], BSobj)



showOneDMR(dmlTest.sm_sig_dmr['1541',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['408',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['2026',], BSobj)
showOneDMR(dmlTest.sm_sig_dmr['6002',], BSobj)


write.table(dmrs.sm2, "dmrs.sm2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
write.table(dmrs.sm, "dmrs.sm.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


write.table(dmlTest.sm_sig_dmr, "dmr_sig.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)





dmr.sm_pmd <- fread("dmrs.sm_pmd.bedGraph")

ggplot(dmr.sm_pmd,aes(x=V16,y=V17,fill=V5))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)

dmr.sm2_pmd <- fread("dmrs.sm2_pmd.bedGraph")

ggplot(dmr.sm2_pmd,aes(x=V16,y=V17,fill=V5))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("HMD","PMD")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )+
  labs(x = "Region Type", 
       y= "Length",
       title="PMD")

dmr.sm_H3K9me3 <- fread("dmrs.sm_H3K9me3.bedGraph")

ggplot(dmr.sm_H3K9me3,aes(x=V20,y=V21,fill=V20))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )+
  labs(x = "Region Type", 
       y= "Length",
       title="PMD")

dmr.sm2_H3K9me3 <- fread("dmrs.sm2_H3K9me3.bedGraph")

ggplot(dmr.sm2_H3K9me3,aes(x=V20,y=V21,fill=V20))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )+
  labs(x = "Region Type", 
       y= "Length",
       title="H3K9me3")+
  geom_jitter(width = 0.05,size=0.5)


dmr.sm2_H3K27me3 <- fread("dmrs.sm2_H3K27me3.bedGraph")

ggplot(dmr.sm2_H3K27me3,aes(x=V20,y=V21,fill=V20))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )+
  labs(x = "Region Type", 
       y= "Length",
       title="H3K27me3")+
  geom_jitter(width = 0.05,size=0.5)


dmr.sm2_H3K36me3 <- fread("dmrs.sm2_H3K36me3.bedGraph")

ggplot(dmr.sm2_H3K36me3,aes(x=V20,y=V21,fill=V20))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )+
  labs(x = "Region Type", 
       y= "Length",
       title="H3K36me3")+
  geom_jitter(width = 0.05,size=0.5)


dmr.sm2_H3K27ac <- fread("dmrs.sm2_H3K27ac.bedGraph")

ggplot(dmr.sm2_H3K27ac,aes(x=V20,y=V21,fill=V20))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )+
  labs(x = "Region Type", 
       y= "Length",
       title="H3K27ac")+
  geom_jitter(width = 0.05,size=0.5)



dmr.sm2_pmd_sub <- subset(dmr.sm2_pmd,subset=V5!="NA")
dmr.sm2_pmd_up <- subset(dmr.sm2_pmd_sub,subset=V16=="Up-regulated")
dmr.sm2_pmd_down <- subset(dmr.sm2_pmd_sub,subset=V16=="Down-regulated")

ggplot(dmr.sm2_pmd_sub,aes(x=V5,y=V17,fill=V16))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 300 )

ggplot(dmr.sm2_pmd_up,aes(x=V5,y=V17,fill=V5))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("HMD","PMD")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )

ggplot(dmr.sm2_pmd_down,aes(x=V5,y=V17,fill=V5))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("HMD","PMD")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 400 )



dmr.sm2_pmd_p <- subset(dmr.sm2_pmd_sub,subset=V5=="PMD")
ggplot(dmr.sm2_pmd_p,aes(x=V16,y=V17,fill=V16))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,800)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 600 )+
  labs(x = "Region Type", 
       y= "Length",
       title="PMD")+
  geom_jitter(width = 0.05,size=0.5)
dmr.sm2_pmd_h <- subset(dmr.sm2_pmd_sub,subset=V5=="HMD")
ggplot(dmr.sm2_pmd_h,aes(x=V16,y=V17,fill=V16))+
  geom_boxplot(outlier.alpha = 0)+ylim(0,600)+
  stat_compare_means(comparisons = list(c("Up-regulated","Down-regulated")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 500 )+
labs(x = "Region Type", 
     y= "Length",
     title="HMD")+
  geom_jitter(width = 0.05,size=0.5)



df_color <- c("#1f76b6", "#ff7d0e")
ggplot(dmrs.sm2.prop, aes(x = 1, y = Freq, fill = Var2))+
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

dmls.sm2.prop <- as.data.frame(t(table(dmls.sm2$state)))
ggplot(dmls.sm2.prop, aes(x = 1, y = Freq, fill = Var2))+
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

h3k9me3_fisherdata <- matrix(c(228, 334, 3799, 5996), nrow = 2)
colnames(h3k9me3_fisherdata) <- c("up", "down") 
rownames(h3k9me3_fisherdata) <- c("k9me3", "other")
fisher.test(h3k9me3_fisherdata)  

h3k27ac_fisherdata <- matrix(c(493, 1352, 3534, 4979), nrow = 2)
colnames(h3k27ac_fisherdata) <- c("up", "down") 
rownames(h3k27ac_fisherdata) <- c("k27ac", "other")
fisher.test(h3k27ac_fisherdata) 




dmr_sig_pmd <- fread("dmr_sig_pmd.bedGraph")
dmr_sig_H3K9me3 <- fread("dmr_sig_H3K9me3.bedGraph")
dmr_sig_H3K27me3 <- fread("dmr_sig_H3K27me3.bedGraph")
dmr_sig_H3K36me3 <- fread("dmr_sig_H3K36me3.bedGraph")
dmr_sig_H3K27ac <- fread("dmr_sig_H3K27ac.bedGraph")
dmr_sig_pmd_p <- subset(dmr_sig_pmd,subset=V5 == 'PMD')
dmr_sig_pmd_h <- subset(dmr_sig_pmd,subset=V5 == 'HMD')
dmr_sig_huh_H3K27ac <- fread("dmr_sig_huh7_H3K27ac.bedGraph")


###total
table(dmlTest.sm_sig_dmr$state)
#Down-regulated   Up-regulated 
#902           2006 

###pmd-hmd
table(dmr_sig_pmd_p$V16)
#Down-regulated   Up-regulated 
#390            783 

table(dmr_sig_pmd_h$V16)
#Down-regulated   Up-regulated 
#464           1150 
sig_pmd_fisherdata <- matrix(c(390, 783, 464, 1150), nrow = 2)
colnames(sig_pmd_fisherdata) <- c("up", "down") 
rownames(sig_pmd_fisherdata) <- c("pmd", "hmd")
fisher.test(sig_pmd_fisherdata)  

sig_pmd_p_fisherdata <- matrix(c(390, 783, 519, 1223), nrow = 2)
colnames(sig_pmd_p_fisherdata) <- c("up", "down") 
rownames(sig_pmd_p_fisherdata) <- c("pmd", "other")
fisher.test(sig_pmd_p_fisherdata)  


sig_pmd_h_fisherdata <- matrix(c(464, 1150, 438, 856), nrow = 2)
colnames(sig_pmd_h_fisherdata) <- c("up", "down") 
rownames(sig_pmd_h_fisherdata) <- c("hmd", "other")
fisher.test(sig_pmd_h_fisherdata)  

###H3K9me3
table(dmr_sig_H3K9me3$V20)
#Down-regulated   Up-regulated 
#25             58
#total
#Down-regulated   Up-regulated 
#902           2006 

sig_H3K9me3_fisherdata <- matrix(c(25, 58, 877, 1948), nrow = 2)
colnames(sig_H3K9me3_fisherdata) <- c("up", "down") 
rownames(sig_H3K9me3_fisherdata) <- c("H3K9me3", "other")
fisher.test(sig_H3K9me3_fisherdata) 




###H3K27me3
table(dmr_sig_H3K27me3$V20)
#Down-regulated   Up-regulated 
#44            78
#total
#Down-regulated   Up-regulated 
#902           2006 

sig_H3K27me3_fisherdata <- matrix(c(44, 78, 858, 1928), nrow = 2)
colnames(sig_H3K27me3_fisherdata) <- c("up", "down") 
rownames(sig_H3K27me3_fisherdata) <- c("H3K27me3", "other")
fisher.test(sig_H3K27me3_fisherdata) 






###H3K27ac
table(dmr_sig_H3K27ac$V14)
#Down-regulated   Up-regulated 
#220           755
#total
#Down-regulated   Up-regulated 
#902           2006 

sig_H3K27ac_fisherdata <- matrix(c(220, 755, 682, 1251), nrow = 2)
colnames(sig_H3K27ac_fisherdata) <- c("up", "down") 
rownames(sig_H3K27ac_fisherdata) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac_fisherdata) 


dmr_sig_H3K27ac$region <- paste(dmr_sig_H3K27ac$V5,dmr_sig_H3K27ac$V6,sep="_")
dmlTest.sm_sig_dmr$region <- paste(dmlTest.sm_sig_dmr$chr,dmlTest.sm_sig_dmr$start,sep = "_")


dmr_sig_huh_H3K27ac0 <- subset(dmr_sig_H3K27ac,subset=V4>5)
dmr_sig_huh_H3K27ac2 <- subset(dmlTest.sm_sig_dmr,subset=region %in% unique(dmr_sig_huh_H3K27ac0$region))
  
table(dmr_sig_huh_H3K27ac2$state)
#Down-regulated   Up-regulated 
# 431      1209
#total
#Down-regulated   Up-regulated 
#902           2006 
sig_H3K27ac2_fisherdata <- matrix(c(431,1209,471,797), nrow = 2)
colnames(sig_H3K27ac2_fisherdata) <- c("up", "down") 
rownames(sig_H3K27ac2_fisherdata) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac2_fisherdata) 

###H3K36me3
table(dmr_sig_H3K36me3$V20)
#Down-regulated   Up-regulated 
#8       18
#total
#Down-regulated   Up-regulated 
#902           2006 

sig_H3K36me3_fisherdata <- matrix(c(8, 18, 894, 1988), nrow = 2)
colnames(sig_H3K36me3_fisherdata) <- c("up", "down") 
rownames(sig_H3K36me3_fisherdata) <- c("H3K36me3", "other")
fisher.test(sig_H3K36me3_fisherdata) 




dmr_sig_H3K4me1 <- fread("dmr_sig_H3K4me1.bedGraph")
dmr_sig_H3K4me1$region <- paste(dmr_sig_H3K4me1$V5,dmr_sig_H3K4me1$V6,sep="_")
dmlTest.sm_sig_dmr$region <- paste(dmlTest.sm_sig_dmr$chr,dmlTest.sm_sig_dmr$start,sep = "_")


dmr_sig_H3K4me1_sub <- subset(dmr_sig_H3K4me1,subset=V4>5)
dmr_sig_H3K4me1_sub2 <- subset(dmlTest.sm_sig_dmr,subset=region %in% unique(dmr_sig_H3K4me1_sub$region))

table(dmr_sig_H3K4me1_sub2$state)
#Down-regulated   Up-regulated 
# 302      617
#total
#Down-regulated   Up-regulated 
#902           2006 
sig_H3K4me1_fisherdata <- matrix(c(302,617,600,1389), nrow = 2)
colnames(sig_H3K4me1_fisherdata) <- c("up", "down") 
rownames(sig_H3K4me1_fisherdata) <- c("H3K4me1", "other")
fisher.test(sig_H3K4me1_fisherdata) 










dmr_sig_H3K9ac <- fread("dmr_sig_H3K9ac.bedGraph")
dmr_sig_H3K9ac$region <- paste(dmr_sig_H3K9ac$V5,dmr_sig_H3K9ac$V6,sep="_")
dmlTest.sm_sig_dmr$region <- paste(dmlTest.sm_sig_dmr$chr,dmlTest.sm_sig_dmr$start,sep = "_")


dmr_sig_H3K9ac_sub <- subset(dmr_sig_H3K9ac,subset=V4>5)
dmr_sig_H3K9ac_sub2 <- subset(dmlTest.sm_sig_dmr,subset=region %in% unique(dmr_sig_H3K9ac_sub$region))

table(dmr_sig_H3K9ac_sub2$state)
#Down-regulated   Up-regulated 
# 502      1352
#total
#Down-regulated   Up-regulated 
#902           2006 
sig_H3K9ac_fisherdata <- matrix(c(502,1352,400,654), nrow = 2)
colnames(sig_H3K9ac_fisherdata) <- c("up", "down") 
rownames(sig_H3K9ac_fisherdata) <- c("H3K9ac", "other")
fisher.test(sig_H3K9ac_fisherdata) 



dmr_sig_H3K36ac <- fread("dmr_sig_H3K36ac.bedGraph")
dmr_sig_H3K36ac$region <- paste(dmr_sig_H3K36ac$V5,dmr_sig_H3K36ac$V6,sep="_")
dmr_sig_H3K36ac_sub <- subset(dmr_sig_H3K36ac,subset=V4>5)
dmr_sig_H3K36ac_sub2 <- subset(dmlTest.sm_sig_dmr,subset=region %in% unique(dmr_sig_H3K36ac_sub$region))

table(dmr_sig_H3K36ac_sub2$state)
#Down-regulated   Up-regulated 
# 449      1229
#total
#Down-regulated   Up-regulated 
#902           2006 
sig_H3K36ac_fisherdata <- matrix(c(449,1229,453,777), nrow = 2)
colnames(sig_H3K36ac_fisherdata) <- c("up", "down") 
rownames(sig_H3K36ac_fisherdata) <- c("H3K36ac", "other")
fisher.test(sig_H3K36ac_fisherdata) 


dmr_sig_H4K5ac <- fread("dmr_sig_H4K5ac.bedGraph")
dmr_sig_H4K5ac$region <- paste(dmr_sig_H4K5ac$V5,dmr_sig_H4K5ac$V6,sep="_")
dmr_sig_H4K5ac_sub <- subset(dmr_sig_H4K5ac,subset=V4>5)
dmr_sig_H4K5ac_sub2 <- subset(dmlTest.sm_sig_dmr,subset=region %in% unique(dmr_sig_H4K5ac_sub$region))

table(dmr_sig_H4K5ac_sub2$state)
#Down-regulated   Up-regulated 
# 331      969
#total
#Down-regulated   Up-regulated 
#902           2006 
sig_H4K5ac_fisherdata <- matrix(c(331,969,571,1037), nrow = 2)
colnames(sig_H4K5ac_fisherdata) <- c("up", "down") 
rownames(sig_H4K5ac_fisherdata) <- c("H4K5ac", "other")
fisher.test(sig_H4K5ac_fisherdata) 








hg38_tss <- fread("hg38_refseq_TSS.bed")
hg38_promoter <- hg38_tss
hg38_promoter$V2 <- hg38_promoter$V2-2000
hg38_promoter$V3 <- hg38_promoter$V3+30
colnames(hg38_promoter) <- c("chr","start","end","length","gene","di")
hg38_promoter <- subset(hg38_promoter,subset=start>0)
write.table(hg38_promoter, "hg38_promoter.bed",sep = "\t",quote=F,row.names = F,col.names = F)



dmr_sig_promoter <- fread("dmr_sig_promoter.bed")
dmr_sig_promoter$region <- paste(dmr_sig_promoter$V7,dmr_sig_promoter$V8,sep="_")



dmr_sig_promoter_sub <- subset(dmlTest.sm_sig_dmr,subset=region %in% unique(dmr_sig_promoter$region))
 

table(dmr_sig_promoter_sub$state)
#Down-regulated   Up-regulated 
# 270      773
#total
#Down-regulated   Up-regulated 
#902           2006 
sig_promoter_fisherdata <- matrix(c(270,773,902,2006), nrow = 2)
colnames(sig_promoter_fisherdata) <- c("up", "down") 
rownames(sig_promoter_fisherdata) <- c("promoter", "other")
fisher.test(sig_promoter_fisherdata) 




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



