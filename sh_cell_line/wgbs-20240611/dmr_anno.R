setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240611/methy_mtx/dmr_anno")
library(data.table)
library(ggrepel)
rm(list = ls())
load("SNG6_KD_Huh7.Rdata")
save.image("SNG6_KD_Huh7.Rdata")


dmr_sig <- fread("SNG6_dmrs.bedGraph")
dmr_sig$V11 <- paste(dmr_sig$V1,dmr_sig$V2,sep = "_")

dmr_pmd <- fread("dmr_pmd.bedGraph")

dmr_pmd_sub <- subset(dmr_pmd,subset=V6!="Neither")

dmr_pmd_sub_p <- subset(dmr_pmd_sub,subset=V6=="commonPMD")

dmr_pmd_sub_h <- subset(dmr_pmd_sub,subset=V6=="commonHMD")


###total
sum_data <- as.data.frame(table(dmr_sig$V10))
#Down-regulated   Up-regulated 
# 49            103

###pmd-hmd
pmd_data <- as.data.frame(table(dmr_pmd_sub_p$V16))
#Down-regulated   Up-regulated 
#  26             40 

hmd_data <- as.data.frame(table(dmr_pmd_sub_h$V16))
#Down-regulated   Up-regulated 
#  13             28 


sig_pmd_p_fisher_data <- as.data.frame(matrix(rbind(pmd_data$Freq,sum_data$Freq),nrow = 2))
sig_pmd_p_fisher_data <- dplyr::select(sig_pmd_p_fisher_data,2,1)
sig_pmd_p_fisher_data[2,] <- sig_pmd_p_fisher_data[2,]-sig_pmd_p_fisher_data[1,]
colnames(sig_pmd_p_fisher_data) <- c("up", "down") 
rownames(sig_pmd_p_fisher_data) <- c("pmd", "other")

fisher.test(sig_pmd_p_fisher_data) 
chisq.test(sig_pmd_p_fisher_data)
sig_pmd_p_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_p_fisher_data)))
sig_pmd_p_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_p_fisher_data)))
sig_pmd_p_result <- matrix(c(sig_pmd_p_chisq_result[1:3,],sig_pmd_p_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_p_result) <- 'PMD'
row.names(sig_pmd_p_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 0.1163 odds ratio  0.5638323 
#X-squared = 2.1869, df = 1, p-value = 0.1392


sig_pmd_h_fisher_data <- as.data.frame(matrix(rbind(hmd_data$Freq,sum_data$Freq),nrow = 2))
sig_pmd_h_fisher_data <- dplyr::select(sig_pmd_h_fisher_data,2,1)
sig_pmd_h_fisher_data[2,] <- sig_pmd_h_fisher_data[2,]-sig_pmd_h_fisher_data[1,]
colnames(sig_pmd_h_fisher_data) <- c("up", "down") 
rownames(sig_pmd_h_fisher_data) <- c("pmd", "other")

fisher.test(sig_pmd_h_fisher_data) 
chisq.test(sig_pmd_h_fisher_data)
sig_pmd_h_chisq_result <- as.data.frame(unlist(chisq.test(sig_pmd_h_fisher_data)))
sig_pmd_h_fisher_result <- as.data.frame(unlist(fisher.test(sig_pmd_h_fisher_data)))
sig_pmd_h_result <- matrix(c(sig_pmd_h_chisq_result[1:3,],sig_pmd_h_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_pmd_h_result) <- 'HMD'
row.names(sig_pmd_h_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value = 1  odds ratio   1.033625
#X-squared = 4.0139e-30, df = 1, p-value = 1



sig_pmd_pvh_fisher_data <- as.data.frame(matrix(rbind(pmd_data$Freq,hmd_data$Freq),nrow = 2))
sig_pmd_pvh_fisher_data <- dplyr::select(sig_pmd_pvh_fisher_data,2,1)

colnames(sig_pmd_pvh_fisher_data) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")
fisher.test(sig_pmd_pvh_fisher_data) 
mcnemar.test(sig_pmd_pvh_fisher_data, correct=F) 

#p-value = 1  odds ratio   0.9981356  
#McNemar's chi-squared = 5.8182, df = 1, p-value = 0.01586












dmr_sig_H4K5ac <- fread("dmr_H4K5ac.bedGraph")
table(dmr_sig_H4K5ac$V14)
dmr_sig_H4K5ac_sub <- subset(dmr_sig_H4K5ac,subset=V4>3)
dmr_sig_H4K5ac_sub$region <- paste(dmr_sig_H4K5ac_sub$V5,dmr_sig_H4K5ac_sub$V6,sep="_")
dmr_sig_H4K5ac_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H4K5ac_sub$region))
H4K5ac_data <- as.data.frame(table(dmr_sig_H4K5ac_sub$V10))
sig_H4K5ac_fisher_data <- as.data.frame(matrix(rbind(H4K5ac_data$Freq,sum_data$Freq),nrow = 2))
sig_H4K5ac_fisher_data <- dplyr::select(sig_H4K5ac_fisher_data,2,1)
sig_H4K5ac_fisher_data[2,] <- sig_H4K5ac_fisher_data[2,]-sig_H4K5ac_fisher_data[1,]
colnames(sig_H4K5ac_fisher_data) <- c("up", "down") 
rownames(sig_H4K5ac_fisher_data) <- c("H4K5ac", "other")
fisher.test(sig_H4K5ac_fisher_data) 
chisq.test(sig_H4K5ac_fisher_data)
sig_H4K5ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H4K5ac_fisher_data)))
sig_H4K5ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H4K5ac_fisher_data)))
sig_H4K5ac_result <- matrix(c(sig_H4K5ac_chisq_result[1:3,],sig_H4K5ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H4K5ac_result) <- 'H4K5ac'
row.names(sig_H4K5ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")
#p-value =  0.5934 odds ratio 0.8837814  
#X-squared = 0.21953, df = 1, p-value = 0.6394







dmr_sig_H3K9ac <- fread("dmr_H3K9ac.bedGraph")
table(dmr_sig_H3K9ac$V14)
dmr_sig_H3K9ac_sub <- subset(dmr_sig_H3K9ac,subset=V4>3)
dmr_sig_H3K9ac_sub$region <- paste(dmr_sig_H3K9ac_sub$V5,dmr_sig_H3K9ac_sub$V6,sep="_")
dmr_sig_H3K9ac_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K9ac_sub$region))
H3K9ac_data <- as.data.frame(table(dmr_sig_H3K9ac_sub$V10))
sig_H3K9ac_fisher_data <- as.data.frame(matrix(rbind(H3K9ac_data$Freq,sum_data$Freq),nrow = 2))
sig_H3K9ac_fisher_data <- dplyr::select(sig_H3K9ac_fisher_data,2,1)
sig_H3K9ac_fisher_data[2,] <- sig_H3K9ac_fisher_data[2,]-sig_H3K9ac_fisher_data[1,]
colnames(sig_H3K9ac_fisher_data) <- c("up", "down") 
rownames(sig_H3K9ac_fisher_data) <- c("H3K9ac", "other")
fisher.test(sig_H3K9ac_fisher_data) 
chisq.test(sig_H3K9ac_fisher_data)
sig_H3K9ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K9ac_fisher_data)))
sig_H3K9ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K9ac_fisher_data)))
sig_H3K9ac_result <- matrix(c(sig_H3K9ac_chisq_result[1:3,],sig_H3K9ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K9ac_result) <- 'H3K9ac'
row.names(sig_H3K9ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")




dmr_sig_H3K27ac <- fread("dmr_H3K27ac.bedGraph")
table(dmr_sig_H3K27ac$V14)
dmr_sig_H3K27ac_sub <- subset(dmr_sig_H3K27ac,subset=V4>3)
dmr_sig_H3K27ac_sub$region <- paste(dmr_sig_H3K27ac_sub$V5,dmr_sig_H3K27ac_sub$V6,sep="_")
dmr_sig_H3K27ac_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K27ac_sub$region))
H3K27ac_data <- as.data.frame(table(dmr_sig_H3K27ac_sub$V10))
sig_H3K27ac_fisher_data <- as.data.frame(matrix(rbind(H3K27ac_data$Freq,sum_data$Freq),nrow = 2))
sig_H3K27ac_fisher_data <- dplyr::select(sig_H3K27ac_fisher_data,2,1)
sig_H3K27ac_fisher_data[2,] <- sig_H3K27ac_fisher_data[2,]-sig_H3K27ac_fisher_data[1,]
colnames(sig_H3K27ac_fisher_data) <- c("up", "down") 
rownames(sig_H3K27ac_fisher_data) <- c("H3K27ac", "other")
fisher.test(sig_H3K27ac_fisher_data) 
chisq.test(sig_H3K27ac_fisher_data)
sig_H3K27ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K27ac_fisher_data)))
sig_H3K27ac_result <- matrix(c(sig_H3K27ac_chisq_result[1:3,],sig_H3K27ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K27ac_result) <- 'H3K27ac'
row.names(sig_H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")





dmr_sig_H3K36ac <- fread("dmr_H3K36ac.bedGraph")
table(dmr_sig_H3K36ac$V14)
dmr_sig_H3K36ac_sub <- subset(dmr_sig_H3K36ac,subset=V4>3)
dmr_sig_H3K36ac_sub$region <- paste(dmr_sig_H3K36ac_sub$V5,dmr_sig_H3K36ac_sub$V6,sep="_")
dmr_sig_H3K36ac_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K36ac_sub$region))
H3K36ac_data <- as.data.frame(table(dmr_sig_H3K36ac_sub$V10))
sig_H3K36ac_fisher_data <- as.data.frame(matrix(rbind(H3K36ac_data$Freq,sum_data$Freq),nrow = 2))
sig_H3K36ac_fisher_data <- dplyr::select(sig_H3K36ac_fisher_data,2,1)
sig_H3K36ac_fisher_data[2,] <- sig_H3K36ac_fisher_data[2,]-sig_H3K36ac_fisher_data[1,]
colnames(sig_H3K36ac_fisher_data) <- c("up", "down") 
rownames(sig_H3K36ac_fisher_data) <- c("H3K36ac", "other")
fisher.test(sig_H3K36ac_fisher_data) 
chisq.test(sig_H3K36ac_fisher_data)
sig_H3K36ac_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K36ac_fisher_data)))
sig_H3K36ac_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K36ac_fisher_data)))
sig_H3K36ac_result <- matrix(c(sig_H3K36ac_chisq_result[1:3,],sig_H3K36ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K36ac_result) <- 'H3K36ac'
row.names(sig_H3K36ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")





dmr_sig_H3K4me1 <- fread("dmr_H3K4me1.bedGraph")
table(dmr_sig_H3K4me1$V14)
dmr_sig_H3K4me1_sub <- subset(dmr_sig_H3K4me1,subset=V4>3)
dmr_sig_H3K4me1_sub$region <- paste(dmr_sig_H3K4me1_sub$V5,dmr_sig_H3K4me1_sub$V6,sep="_")
dmr_sig_H3K4me1_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_H3K4me1_sub$region))
H3K4me1_data <- as.data.frame(table(dmr_sig_H3K4me1_sub$V10))
sig_H3K4me1_fisher_data <- as.data.frame(matrix(rbind(H3K4me1_data$Freq,sum_data$Freq),nrow = 2))
sig_H3K4me1_fisher_data <- dplyr::select(sig_H3K4me1_fisher_data,2,1)
sig_H3K4me1_fisher_data[2,] <- sig_H3K4me1_fisher_data[2,]-sig_H3K4me1_fisher_data[1,]
colnames(sig_H3K4me1_fisher_data) <- c("up", "down") 
rownames(sig_H3K4me1_fisher_data) <- c("H3K4me1", "other")
fisher.test(sig_H3K4me1_fisher_data) 
chisq.test(sig_H3K4me1_fisher_data)
sig_H3K4me1_chisq_result <- as.data.frame(unlist(chisq.test(sig_H3K4me1_fisher_data)))
sig_H3K4me1_fisher_result <- as.data.frame(unlist(fisher.test(sig_H3K4me1_fisher_data)))
sig_H3K4me1_result <- matrix(c(sig_H3K4me1_chisq_result[1:3,],sig_H3K4me1_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_H3K4me1_result) <- 'H3K4me1'
row.names(sig_H3K4me1_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")













###promoter
dmr_sig_promoter <- fread("dmr_promoter.bedGraph")
dmr_sig_promoter$region <- paste(dmr_sig_promoter$V7,dmr_sig_promoter$V8,sep="_")
dmr_sig_promoter_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_promoter$region))
dmr_promoter_data <- as.data.frame(table(dmr_sig_promoter_sub$V10))
dmr_sig_promoter_fisher_data <- as.data.frame(matrix(rbind(dmr_promoter_data$Freq,sum_data$Freq),nrow = 2))
dmr_sig_promoter_fisher_data <- dplyr::select(dmr_sig_promoter_fisher_data,2,1)
dmr_sig_promoter_fisher_data[2,] <- dmr_sig_promoter_fisher_data[2,]-dmr_sig_promoter_fisher_data[1,]
colnames(dmr_sig_promoter_fisher_data) <- c("up", "down") 
rownames(dmr_sig_promoter_fisher_data) <- c("promoter", "other")

fisher.test(dmr_sig_promoter_fisher_data)
chisq.test(dmr_sig_promoter_fisher_data)
sig_promoter_chisq_result <- as.data.frame(unlist(chisq.test(dmr_sig_promoter_fisher_data)))
sig_promoter_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_promoter_fisher_data)))
sig_promoter_result <- matrix(c(sig_promoter_chisq_result[1:3,],sig_promoter_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_promoter_result) <- 'promoter'
row.names(sig_promoter_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")


###promoter
dmr_sig_enhancer <- fread("dmr_enhancer.bedGraph")
dmr_sig_enhancer$region <- paste(dmr_sig_enhancer$V7,dmr_sig_enhancer$V8,sep="_")
dmr_sig_enhancer_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_sig_enhancer$region))
dmr_enhancer_data <- as.data.frame(table(dmr_sig_enhancer_sub$V10))
dmr_sig_enhancer_fisher_data <- as.data.frame(matrix(rbind(dmr_enhancer_data$Freq,sum_data$Freq),nrow = 2))
dmr_sig_enhancer_fisher_data <- dplyr::select(dmr_sig_enhancer_fisher_data,2,1)
dmr_sig_enhancer_fisher_data[2,] <- dmr_sig_enhancer_fisher_data[2,]-dmr_sig_enhancer_fisher_data[1,]
colnames(dmr_sig_enhancer_fisher_data) <- c("up", "down") 
rownames(dmr_sig_enhancer_fisher_data) <- c("enhancer", "other")

fisher.test(dmr_sig_enhancer_fisher_data)
chisq.test(dmr_sig_enhancer_fisher_data)
sig_enhancer_chisq_result <- as.data.frame(unlist(chisq.test(dmr_sig_enhancer_fisher_data)))
sig_enhancer_fisher_result <- as.data.frame(unlist(fisher.test(dmr_sig_enhancer_fisher_data)))
sig_enhancer_result <- matrix(c(sig_enhancer_chisq_result[1:3,],sig_enhancer_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(sig_enhancer_result) <- 'enhancer'
row.names(sig_enhancer_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")


dmr_sig_test_result <- as.data.frame(t(
  cbind(sig_promoter_result,sig_enhancer_result,
        sig_H3K4me1_result,sig_H3K27ac_result,
        sig_H3K9ac_result,sig_H4K5ac_result,
        sig_H3K36ac_result,
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

dmr_sig_test_result$case = factor(dmr_sig_test_result$case, levels=c("Enriched in Hypo","Enriched in Hyper"))
ggplot(dmr_sig_test_result,aes(x=lnR,y=logfdr))+
  geom_point(size = 1,aes(color = case))+# 根据expression水平进行着色
  xlab(expression("ln odds ratio")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  theme_bw() +
  geom_text_repel(aes(label=histone_mark, color = case), size =5,hjust=1,vjust=-1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Chi-square test result of Histone Mark in SNHG6-KD-Huh7")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")



dmr_sig_bed_down <- subset(dmr_sig,subset=V9 >0)
dmr_sig_bed_up <- subset(dmr_sig,subset=V9 <0)
dmr_sig_bed_down <- dplyr::select(dmr_sig_bed_down,1,2,3)
dmr_sig_bed_down$V4 <- row.names(dmr_sig_bed_down)

dmr_sig_bed_up <- dplyr::select(dmr_sig_bed_up,1,2,3)
dmr_sig_bed_up$V4 <- row.names(dmr_sig_bed_up)

write.table(dmr_sig_bed_up,"dmr_up.bed",sep = "\t",quote=F,row.names = F,col.names = F)
write.table(dmr_sig_bed_down,"dmr_down.bed",sep = "\t",quote=F,row.names = F,col.names = F)

up_stat <- fread("up_stat.txt")
up_anno <- fread("up_annotation.txt")
up_stat <- up_stat[15:46,] 
up_stat$dmr_type <- 'up'
down_stat <- fread("down_stat.txt")
down_anno <- fread("down_annotation.txt")
down_stat <- down_stat[14:46,] 
down_stat$dmr_type <- 'down'


stat_sum <- rbind(up_stat,down_stat)
colnames(stat_sum) <- c('Annotation','Nums','TotalSize','Log2Ratio','LogP','dmr_type')
stat_sum$Nums <- as.numeric(stat_sum$Nums)
stat_sum$Log2Ratio <- as.numeric(stat_sum$Log2Ratio)
stat_sum <- subset(stat_sum,subset=Nums!=0)

ggplot(stat_sum,aes(x=Annotation,y=dmr_type))+
  geom_point(aes(size=Nums,color=Log2Ratio))+
  scale_color_gradient(low = "#73BEF5",high = "#F2055C")+theme_bw()+
  theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),legend.text=element_text(size = 10))+
  theme(axis.text.x = element_text(size = 14,color="black",angle=45,hjust = 1),axis.text.y = element_text(size = 15,color="black"))+
  theme(panel.grid = element_blank())

write.table(stat_sum,"stat_sum.txt")