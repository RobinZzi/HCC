setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240705/bismark_cov/dmr_dmcp0.01delta0.2")
library(Rcircos)

install.packages("RCircos")
rm(list = ls())
save.image("GA45_KD_293t.Rdata")
dmr_sig <- fread("GA45_dmrs.bedGraph")
dmr_sig$V11 <- paste(dmr_sig$V1,dmr_sig$V2,sep = "_")

dmr_pmd <- fread("dmr_pmd.bedGraph")

dmr_pmd_sub <- subset(dmr_pmd,subset=V6!="Neither")

dmr_pmd_sub_p <- subset(dmr_pmd_sub,subset=V6=="commonPMD")

dmr_pmd_sub_h <- subset(dmr_pmd_sub,subset=V6=="commonHMD")


###total
sum_data <- as.data.frame(table(dmr_sig$V10))


###pmd-hmd
pmd_data <- as.data.frame(table(dmr_pmd_sub_p$V16))
hmd_data <- as.data.frame(table(dmr_pmd_sub_h$V16))


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




sig_pmd_pvh_fisher_data <- as.data.frame(matrix(rbind(pmd_data$Freq,hmd_data$Freq),nrow = 2))
sig_pmd_pvh_fisher_data <- dplyr::select(sig_pmd_pvh_fisher_data,2,1)

colnames(sig_pmd_pvh_fisher_data) <- c("up", "down") 
rownames(sig_pmd_pvh_fisher_data) <- c("pmd", "hmd")
sig_pmd_pvh_fisher_data <- as.matrix(sig_pmd_pvh_fisher_data)

mcnemar.test(sig_pmd_pvh_fisher_data, correct=F) 






dmr_H3K27ac <- fread("dmr_H3K27ac.bedGraph")
table(dmr_H3K27ac$V14)
dmr_H3K27ac_sub <- subset(dmr_H3K27ac,subset=V4>3)
dmr_H3K27ac_sub$region <- paste(dmr_H3K27ac_sub$V5,dmr_H3K27ac_sub$V6,sep="_")
dmr_H3K27ac_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_H3K27ac_sub$region))
H3K27ac_data <- as.data.frame(table(dmr_H3K27ac_sub$V10))
H3K27ac_fisher_data <- as.data.frame(matrix(rbind(H3K27ac_data$Freq,sum_data$Freq),nrow = 2))


H3K27ac_fisher_data <- dplyr::select(H3K27ac_fisher_data,2,1)
H3K27ac_fisher_data[2,] <- H3K27ac_fisher_data[2,]-H3K27ac_fisher_data[1,]
colnames(H3K27ac_fisher_data) <- c("up", "down") 
rownames(H3K27ac_fisher_data) <- c("H3K27ac", "other")
fisher.test(H3K27ac_fisher_data) 
chisq.test(H3K27ac_fisher_data,correct = T)
H3K27ac_chisq_result <- as.data.frame(unlist(chisq.test(H3K27ac_fisher_data)))
H3K27ac_fisher_result <- as.data.frame(unlist(fisher.test(H3K27ac_fisher_data)))
H3K27ac_result <- matrix(c(H3K27ac_chisq_result[1:3,],H3K27ac_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(H3K27ac_result) <- 'H3K27ac'
row.names(H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")



dmr_H3K4me3 <- fread("dmr_H3K4me3.bedGraph")
table(dmr_H3K4me3$V14)
dmr_H3K4me3_sub <- subset(dmr_H3K4me3,subset=V4>3)
dmr_H3K4me3_sub$region <- paste(dmr_H3K4me3_sub$V5,dmr_H3K4me3_sub$V6,sep="_")
dmr_H3K4me3_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_H3K4me3_sub$region))
H3K4me3_data <- as.data.frame(table(dmr_H3K4me3_sub$V10))
H3K4me3_fisher_data <- as.data.frame(matrix(rbind(H3K4me3_data$Freq,sum_data$Freq),nrow = 2))
H3K4me3_fisher_data <- dplyr::select(H3K4me3_fisher_data,2,1)
H3K4me3_fisher_data[2,] <- H3K4me3_fisher_data[2,]-H3K4me3_fisher_data[1,]
colnames(H3K4me3_fisher_data) <- c("up", "down") 
rownames(H3K4me3_fisher_data) <- c("H3K4me3", "other")
fisher.test(H3K4me3_fisher_data) 
chisq.test(H3K4me3_fisher_data,correct = T)
H3K4me3_chisq_result <- as.data.frame(unlist(chisq.test(H3K4me3_fisher_data)))
H3K4me3_fisher_result <- as.data.frame(unlist(fisher.test(H3K4me3_fisher_data)))
H3K4me3_result <- matrix(c(H3K4me3_chisq_result[1:3,],H3K4me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(H3K4me3_result) <- 'H3K4me3'
row.names(H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")



dmr_H3K4me3_2 <- fread("dmr_H3K4me3_2.bedGraph")
table(dmr_H3K4me3_2$V14)
dmr_H3K4me3_2_sub <- subset(dmr_H3K4me3_2,subset=V4>3)
dmr_H3K4me3_2_sub$region <- paste(dmr_H3K4me3_2_sub$V5,dmr_H3K4me3_2_sub$V6,sep="_")
dmr_H3K4me3_2_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_H3K4me3_2_sub$region))
H3K4me3_2_data <- as.data.frame(table(dmr_H3K4me3_2_sub$V10))
H3K4me3_2_fisher_data <- as.data.frame(matrix(rbind(H3K4me3_2_data$Freq,sum_data$Freq),nrow = 2))
H3K4me3_2_fisher_data <- dplyr::select(H3K4me3_2_fisher_data,2,1)
H3K4me3_2_fisher_data[2,] <- H3K4me3_2_fisher_data[2,]-H3K4me3_2_fisher_data[1,]
colnames(H3K4me3_2_fisher_data) <- c("up", "down") 
rownames(H3K4me3_2_fisher_data) <- c("H3K4me3_2", "other")
fisher.test(H3K4me3_2_fisher_data) 
chisq.test(H3K4me3_2_fisher_data,correct = T)
H3K4me3_2_chisq_result <- as.data.frame(unlist(chisq.test(H3K4me3_2_fisher_data)))
H3K4me3_2_fisher_result <- as.data.frame(unlist(fisher.test(H3K4me3_2_fisher_data)))
H3K4me3_2_result <- matrix(c(H3K4me3_2_chisq_result[1:3,],H3K4me3_2_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(H3K4me3_2_result) <- 'H3K4me3_2'
row.names(H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")




dmr_promoter <- fread("dmr_promoter.bedGraph")
promoter_data <- as.data.frame(table(dmr_promoter$V16))
promoter_fisher_data <- as.data.frame(matrix(rbind(promoter_data$Freq,sum_data$Freq),nrow = 2))
promoter_fisher_data <- dplyr::select(promoter_fisher_data,2,1)
promoter_fisher_data[2,] <- promoter_fisher_data[2,]-promoter_fisher_data[1,]
colnames(promoter_fisher_data) <- c("up", "down") 
rownames(promoter_fisher_data) <- c("promoter", "other")
fisher.test(promoter_fisher_data) 
chisq.test(promoter_fisher_data,correct = T)
promoter_chisq_result <- as.data.frame(unlist(chisq.test(promoter_fisher_data)))
promoter_fisher_result <- as.data.frame(unlist(fisher.test(promoter_fisher_data)))
promoter_result <- matrix(c(promoter_chisq_result[1:3,],promoter_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(promoter_result) <- 'promoter'
row.names(promoter_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")


dmr_H3K9me3 <- fread("dmr_H3K9me3.bedGraph")
table(dmr_H3K9me3$V14)
dmr_H3K9me3_sub <- subset(dmr_H3K9me3,subset=V4>3)
dmr_H3K9me3_sub$region <- paste(dmr_H3K9me3_sub$V5,dmr_H3K9me3_sub$V6,sep="_")
dmr_H3K9me3_sub <- subset(dmr_sig,subset=V11 %in% unique(dmr_H3K9me3_sub$region))
H3K9me3_data <- as.data.frame(table(dmr_H3K9me3_sub$V10))
H3K9me3_fisher_data <- as.data.frame(matrix(rbind(H3K9me3_data$Freq,sum_data$Freq),nrow = 2))
H3K9me3_fisher_data <- dplyr::select(H3K9me3_fisher_data,2,1)
H3K9me3_fisher_data[2,] <- H3K9me3_fisher_data[2,]-H3K9me3_fisher_data[1,]
colnames(H3K9me3_fisher_data) <- c("up", "down") 
rownames(H3K9me3_fisher_data) <- c("H3K9me3", "other")
fisher.test(H3K9me3_fisher_data) 
chisq.test(H3K9me3_fisher_data,correct = T)
H3K9me3_chisq_result <- as.data.frame(unlist(chisq.test(H3K9me3_fisher_data)))
H3K9me3_fisher_result <- as.data.frame(unlist(fisher.test(H3K9me3_fisher_data)))
H3K9me3_result <- matrix(c(H3K9me3_chisq_result[1:3,],H3K9me3_fisher_result[c('p.value','estimate.odds ratio'),]),nrow = 5)
colnames(H3K9me3_result) <- 'H3K9me3'
row.names(H3K27ac_result) <- c("X_squared","df","chisq_p","fisher_p","odds_ratio")





dmr_sig_test_result <- as.data.frame(t(
  cbind(sig_pmd_h_result,sig_pmd_p_result,promoter_result,
        H3K27ac_result,
        H3K9me3_result,
        H3K4me3_result,H3K4me3_2_result)))
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
  geom_text_repel(aes(label=histone_mark, color = case), size =5,hjust=1,vjust=1)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Chi-square test result of Histone Mark in GADD45A-KD-293T")+
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




up_stat <- fread("stat.txt")
up_anno <- fread("annotation.txt")
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
stat_sum$LogP <- -stat_sum$LogP 
stat_sum <- subset(stat_sum,subset=Nums!=0)

ggplot(stat_sum,aes(x=Annotation,y=dmr_type))+
  geom_point(aes(size=Nums,color=Log2Ratio))+
  scale_color_gradient(low = "#73BEF5",high = "#F2055C")+theme_bw()+
  theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),legend.text=element_text(size = 10))+
  theme(axis.text.x = element_text(size = 14,color="black",angle=45,hjust = 1),axis.text.y = element_text(size = 15,color="black"))+
  theme(panel.grid = element_blank())

write.table(stat_sum,"stat_sum.txt")
  
