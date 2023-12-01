setwd("/Users/huyiping/desktop/TCGA-LIHC")
setwd("/storage/zhangyanxiaoLab/zhangliwen/Projects/sc-LIHC/scripts&cache_data/TCGA-LIHC")
if(T){
  library(data.table)
  library(pheatmap)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(limma)
  library(DESeq2)
  library(ggstatsplot)
  library(ggside)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(survminer) 
  library(ggcorrplot)
  library(survival)
}
rm(list=ls())

load("TCGA-LIHC.Rdata")
save.image("TCGA-LIHC.Rdata")

meth <- fread("TCGA-LIHC.methylation450.tsv",header = T)
pheno <- fread("TCGA-LIHC.phenotype.tsv",header = T)
anno <- fread("illuminaMethyl450_hg38_GDC")
fpkm <- fread("TCGA-LIHC.htseq_fpkm.tsv")
surv <- fread("TCGA-LIHC.survival.tsv")
count <- fread("TCGA-LIHC.htseq_counts.tsv")


#meth_filt#
meth_filt <-subset(meth,meth$`Composite Element REF` %in% anno$`#id`)
rownames(meth_filt) <- anno$`#id`
meth_filt2 <- na.omit(meth_filt)
meth_filt2 <- meth_filt2[,-1]







#id-transform#
fpkm$Ensembl_ID2 <- fpkm$Ensembl_ID
fpkm$Ensembl_ID2 <- str_sub(fpkm$Ensembl_ID,1,15)
Ensembl_ID <- fpkm$Ensembl_ID2
Ensembl_ID <- str_sub(Ensembl_ID,1,15)
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
fpkm_filt <- subset(fpkm, fpkm$Ensembl_ID2 %in% gene_symbol$ENSEMBL)
gene_symbol= gene_symbol[match(fpkm_filt$Ensembl_ID2,gene_symbol$ENSEMBL),]
fpkm_filt$SYMBOL <- gene_symbol$SYMBOL
exp <- fpkm_filt
exp2<- avereps(exp,ID=exp$SYMBOL)
exp2 <- as.data.frame(exp2)
exp2 <- exp2[!is.na(exp2$SYMBOL),]#去除没有匹配的
exp2 <- exp2[,-1]
row.names(exp2) <- exp2$SYMBOL
exp2 <- exp2[,1:424]


#cobine-exp&meth#
meth_id <-colnames(meth_filt2)
exp_id <- colnames(exp2)
sampleid <- intersect(exp_id,meth_id)
exp3 <- exp2[,sampleid]
meth3 <- as.data.frame(meth_filt2)
meth3 <- meth3[,sampleid]
rownames(meth3) <- meth_filt2$`Composite Element REF`
meth4 <- colMeans(meth3)
meth4 <- as.data.frame(as.numeric(meth4))
colnames(meth4) <- "mean_meth"


exp4 <- as.data.frame(t(exp3))
exp4=as.data.frame(lapply(exp4,as.numeric))

comb <- cbind(meth4,exp4)

row.names(comb) <- sampleid

ggplot(data=comb,aes(x=SNHG6,y=mean_meth))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=comb, method = "pearson")

ggplot(data=comb,aes(x=LSM8,y=mean_meth))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=comb, method = "pearson")

ggplot(data=comb,aes(x=XRN2,y=mean_meth))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=comb, method = "pearson")

ggscatterstats(data = comb, 
               y = SNHG6, 
               x = mean_meth,
               centrality.para = "mean",
)

ggscatterstats(data = comb, 
               y = LSM8,
               x = mean_meth,
               centrality.para = "mean"
)

ggscatterstats(data = comb, 
               y = GADD45A,
               x = mean_meth,
               centrality.para = "mean"
)


ggscatterstats(data = comb, 
               y = GADD45A,
               x = mean_meth,
               centrality.para = "mean"
)





VlnPlot(hpc,features = "UHRF1",group.by = "sample")



tumorid <- grep(pattern = '-01',x=sampleid,value = T)
normalid <- grep(pattern = '-11',x=sampleid,value = T)


row.names(comb) <- sampleid
tumorcomb <- comb[tumorid,]
normalcomb <- comb[normalid,]

tumorcom_log <- log(tumorcomb+1)
y <- as.numeric(tumorcomb[,"mean_meth"])

colnames <- colnames(tumorcomb)

cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  
  test <- cor.test(as.numeric(tumorcomb[,i]),y,type="spearman")
  
  cor_data_df[i,2] <- test$estimate
  
  cor_data_df[i,3] <- test$p.value
  
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
cor_data_df$fdr <- p.adjust(cor_data_df$pvalue,method="fdr",n=length(cor_data_df$pvalue))


write.table(cor_data_df, "TCGA-example-meth_mean-correlation-fdr.txt",sep = "\t",quote=F,row.names = T,col.names = F)


cor_data_sig <- cor_data_df %>% 
  
  filter(fdr < 0.05) %>% 
  
  arrange(desc(abs(correlation)))%>% 
  
  dplyr::slice(1:500)


write.table(cor_data_df_log,"TCGA-example-log-meth_mean-correlation.txt",quote=F,row.names = T,col.names = F)

write.table(tumorcomb,"TCGA-example-exp.txt",quote=F,row.names = T,col.names = F)

ggscatterstats(data = tumorcomb, 
               y = ETS1, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = XRN2, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = VIM, 
               x = mean_meth,
               centrality.para = "mean"
)
ggscatterstats(data = tumorcomb, 
               y = UHRF1, 
               x = mean_meth,
               centrality.para = "mean",
)

ggscatterstats(data = tumorcomb, 
               y = TIMP1, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = mean_meth, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = PRKX, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = SSX3, 
               x = mean_meth,
               centrality.para = "mean",
)

ggscatterstats(data = tumorcomb, 
               y = LDHB, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = SPINT2, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = METTL3, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = FXR1, 
               x = mean_meth,
               centrality.para = "mean",
)
ggscatterstats(data = tumorcomb, 
               y = GADD45A, 
               x = mean_meth,
               centrality.para = "mean",
)




ggscatterstats(data = tumorcomb, 
               y = MET, 
               x = mean_meth,
               centrality.para = "mean",
)




ggscatterstats(data = tumorcomb, 
               y = ONECUT2, 
               x = mean_meth,
               centrality.para = "mean",
)







anno2 <- anno[,-2] 
anno2 <- anno2[,-5]
bigmeth <- cbind(anno2,meth_filt)
bigmeth[,5] <- round(bigmeth[,3]/100000)
bigmeth2 <- bigmeth[,-3]
bigmeth2 <- bigmeth2[,-3]
bigmeth2 <- tidyr::unite(bigmeth2,"id",chrom,'Composite Element REF',sep="_",remove=TRUE)
bigmeth2 <- aggregate(bigmeth2[,3:432],by=list(id=bigmeth2$id),FUN=mean)
bigmeth3 <- bigmeth2[-1,]
bigmeth4 <- na.omit(bigmeth3)
row.names(bigmeth4) <- bigmeth4[,1]
bigmeth4 <- bigmeth4[,-1]



TCGA_pheatmap <- pheatmap(bigmeth4,
                          show_rownames = F, show_colnames = F,
                          cluster_rows = T, cluster_cols = T,
                          clustering_method = "average",
                          color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                          treeheight_row = 100,
                          treeheight_col = 20,
                          fontsize          = 12)

pheatmap(bigmeth4,
         show_rownames = F, show_colnames = F,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 100,
         treeheight_col = 20,
         fontsize          = 12)

pheatmap(bigmeth4[,colanno_row],
         show_rownames = F, show_colnames = F,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 100,
         treeheight_col = 20,
         annotation_col = colanno,
         annotation_colors = ann_colors,
         fontsize          = 12)

TCGA_order_row = TCGA_pheatmap$tree_row$order
TCGA_order_col = TCGA_pheatmap$tree_col$order
sort_TAGA_data = data.frame(bigmeth4[TCGA_order_row,TCGA_order_col])

hypogroup <- colnames(sort_TAGA_data[,1:107])
hypergroup <- colnames(sort_TAGA_data[,112:270])
hypogroup2 <- colnames(sort_TAGA_data[,270:430])
example <- as.data.frame(sort_TAGA_data[,1])

example$cite <- row.names(sort_TAGA_data)
example <- as.data.frame(str_split_fixed(example$cite,"_",2))


colnames(example) <- c("chr","start")
example$start <- as.numeric(example$start)
example$end <- example$start+1
example$start <- example$start*100000
example$end <- example$end*100000
example$end <- as.numeric(example$end)

example$count <- as.numeric(sort_TAGA_data[,1])
example$count <- as.numeric(example$count)
write.table(example, "TCGA-example.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

example2 <- example
example2$count <- as.numeric(sort_TAGA_data[,200])
write.table(example2, "TCGA-example2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


example3 <- example
example3$count <- as.numeric(sort_TAGA_data[,400])
write.table(example3, "TCGA-example3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


bigmeth5  <- lapply(bigmeth3, function(z) {
  mtx <- cbind(c(NA, head(z, -1)), z, c(tail(z, -1), NA))
  mtx[is.na(mtx[,2]) & rowSums(is.na(mtx)) > 1,] <- NA
  out <- ifelse(is.na(mtx[,2]), rowMeans(mtx, na.rm = TRUE), mtx[,2])
  out[is.nan(out)] <- NA
  out
})
bigmeth5 <- as.data.frame(bigmeth5)
bigmeth6 <- na.omit(bigmeth5)
row.names(bigmeth6) <- bigmeth6[,1]
cites2 <- row.names(bigmeth6)
bigmeth6 <- bigmeth6[,-1]
bigmeth6 <- as.data.frame(lapply(bigmeth6,as.numeric))
row.names(bigmeth6) <- cites2
TCGA_pheatmap2 <- pheatmap(bigmeth6,
                           show_rownames = F, show_colnames = F,
                           cluster_rows = T, cluster_cols = T,
                           clustering_method = "average",
                           color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                           treeheight_row = 100,
                           treeheight_col = 20,
                           fontsize          = 12)

bigmeth7 <- bigmeth6

colnames(bigmeth7)<-gsub("\\.","-",colnames(bigmeth7))
tcgaheatmap3 <- pheatmap(bigmeth7[,colanno_row],
                         show_rownames = F, show_colnames = F,
                         cluster_rows = T, cluster_cols = T,
                         clustering_method = "average",
                         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                         cellwidth = 1,cellheight = 0.02,
                         treeheight_row = 100,
                         treeheight_col = 20,
                         annotation_col = colanno,
                         annotation_colors = ann_colors,
                         fontsize          = 12)

pheatmap(bigmeth7[,colanno_row],
         show_rownames = F, show_colnames = F,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         cellwidth = 1,cellheight = 0.02,
         treeheight_row = 100,
         treeheight_col = 20,
         annotation_col = colanno,
         annotation_colors = ann_colors,
         fontsize          = 12)

colanno_tumor <- as.data.frame(tumorid)
colnames(colanno_tumor) <- "id"
colanno_tumor$sample <- "tumor"

colanno_normal <- as.data.frame(normalid)
colnames(colanno_normal) <- "id"
colanno_normal$sample <- "normal"

colanno <- rbind(colanno_tumor,colanno_normal)
rownames(colanno) <- colanno$id
colanno_row <- rownames(colanno)
colanno <- as.data.frame(colanno[,-1])
rownames(colanno) <- colanno_row
colnames(colanno) <- "sample"
ann_colors = list(
  sample = c(normal = "#fff4e1",tumor = "#80d6ff")
)




TCGA_order_row2 = TCGA_pheatmap2$tree_row$order
TCGA_order_col2 = TCGA_pheatmap2$tree_col$order
sort_TAGA_data2 = data.frame(bigmeth6[TCGA_order_row2,TCGA_order_col2])


re <- as.data.frame(str_split_fixed(row.names(sort_TAGA_data2[8000:15996,]),"_",2))
re2 <- as.data.frame(str_split_fixed(row.names(sort_TAGA_data2[11358:15996,]),"_",2))
write.table(re,"re.bed",sep = "\t",quote=F,row.names = F,col.names = F)


inter <- read.table("TCGA.bed",sep = "\t")
combines <- read.table("combine.bed",sep = "\t")
example0 <- as.data.frame(str_split_fixed(row.names(sort_TAGA_data2),"_",2))


colnames(re) <- c("chr","start")
re$start <- as.numeric(re$start)
re$end <- re$start+1
re$start <- re$start*100000
re$end <- re$end*100000
re$end <- as.numeric(re$end)
write.table(re,"re.bed",sep = "\t",quote=F,row.names = F,col.names = F)
colnames(re2) <- c("chr","start")
re2$start <- as.numeric(re2$start)
re2$end <- re2$start+1
re2$start <- re2$start*100000
re2$end <- re2$end*100000
re2$end <- as.numeric(re2$end)
write.table(re2,"re2.bed",sep = "\t",quote=F,row.names = F,col.names = F)


colnames(example0) <- c("chr","start")
example0$start <- as.numeric(example0$start)
example0$end <- example0$start+1
example0$start <- example0$start*100000
example0$end <- example0$end*100000
example0$end <- as.numeric(example0$end)


example4 <- example0
example4$count <- as.numeric(sort_TAGA_data2[,1])
example4$count <- as.numeric(example4$count)
write.table(example4, "TCGA-example4.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

example5 <- example0
example5$count <- as.numeric(sort_TAGA_data2[,200])
write.table(example5, "TCGA-example5.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


example6 <- example0
example6$count <- as.numeric(sort_TAGA_data2[,400])
write.table(example6, "TCGA-example6.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)




cpgs <- anno2[,-1]
write.table(cpgs, "450klocate.bed",sep = "\t",quote=F,row.names = F,col.names = F)

















bigmeth_10k <- cbind(anno2,meth_filt)
bigmeth_10k[,5] <- round(bigmeth_10k[,3]/10000)
bigmeth_10k_2 <- bigmeth_10k[,-3]
bigmeth_10k_2 <- bigmeth_10k_2[,-3]
bigmeth_10k_2 <- tidyr::unite(bigmeth_10k_2,"id",chrom,'Composite Element REF',sep="_",remove=TRUE)
bigmeth_10k_2 <- aggregate(bigmeth_10k_2[,3:432],by=list(id=bigmeth_10k_2$id),FUN=mean)
bigmeth_10k_3 <- bigmeth_10k_2[-1,]
bigmeth_10k_4  <- lapply(bigmeth_10k_3, function(z) {
  mtx <- cbind(c(NA, head(z, -1)), z, c(tail(z, -1), NA))
  mtx[is.na(mtx[,2]) & rowSums(is.na(mtx)) > 1,] <- NA
  out <- ifelse(is.na(mtx[,2]), rowMeans(mtx, na.rm = TRUE), mtx[,2])
  out[is.nan(out)] <- NA
  out
})
bigmeth_10k_4 <- as.data.frame(bigmeth_10k_4)
bigmeth_10k_5 <- na.omit(bigmeth_10k_4)
row.names(bigmeth_10k_5) <- bigmeth_10k_5[,1]
cites3 <- row.names(bigmeth_10k_5)
bigmeth_10k_5 <- bigmeth_10k_5[,-1]
bigmeth_10k_5 <- as.data.frame(lapply(bigmeth_10k_5,as.numeric))
row.names(bigmeth_10k_5) <- cites3
TCGA_pheatmap3 <- pheatmap(bigmeth_10k_5,
                           show_rownames = F, show_colnames = F,
                           cluster_rows = T, cluster_cols = T,
                           clustering_method = "average",
                           color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                           treeheight_row = 100,
                           treeheight_col = 20,
                           fontsize          = 12)


colanno <- as.data.frame()








TCGA_order_row3 = TCGA_pheatmap3$tree_row$order
TCGA_order_col3 = TCGA_pheatmap3$tree_col$order
sort_TAGA_data3 = data.frame(bigmeth6[TCGA_order_row2,TCGA_order_col2])


example0 <- as.data.frame(str_split_fixed(row.names(sort_TAGA_data2),"_",2))


colnames(example0) <- c("chr","start")
example0$start <- as.numeric(example0$start)
example0$end <- example0$start+1
example0$start <- example0$start*100000
example0$end <- example0$end*100000
example0$end <- as.numeric(example0$end)


example4 <- example0
example4$count <- as.numeric(sort_TAGA_data2[,1])
example4$count <- as.numeric(example4$count)
write.table(example4, "TCGA-example4.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

example5 <- example0
example5$count <- as.numeric(sort_TAGA_data2[,200])
write.table(example5, "TCGA-example5.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


example6 <- example0
example6$count <- as.numeric(sort_TAGA_data2[,400])
write.table(example6, "TCGA-example6.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)




cpgs <- anno2[,-1]
write.table(cpgs, "450klocate.bed",sep = "\t",quote=F,row.names = F,col.names = F)


hypermeth_group <- as.data.frame(hypergroup)
hypometh_group2 <- as.data.frame(hypogroup2)
hypometh_group2$level <- "hypo2"
hypermeth_group$level <- "hyper"
colnames(hypermeth_group) <- c("sample","level")
colnames(hypometh_group2) <- c("sample","level")
hypometh_group <- as.data.frame(hypogroup)
hypometh_group$level <- "hypo"
colnames(hypometh_group) <- c("sample","level")
methgroup <- rbind(hypometh_group,hypermeth_group,hypometh_group2)


surv$sample <- gsub("\\-",".",surv$sample)
surv_filt <- subset(surv,surv$sample %in% methgroup$sample)
surv_filt_data <- merge(surv_filt,methgroup,all=TRUE, sort=TRUE)
surv_filt_data <- subset(surv_filt_data,surv_filt_data$sample %in% surv_filt$sample)


surv_filt_data$month <- round(surv_filt_data$OS.time/30,2)
survana <- Surv(time = surv_filt_data$month,event =surv_filt_data$OS=='1')

KMfit <- survfit(survana  ~ surv_filt_data$level)

ggsurvplot(KMfit,                     #拟合对象
           data = surv_filt_data,               #变量数据来源
           pval = TRUE,               #P值
           surv.median.line = "hv",   #中位生存时间线
           risk.table = F,         #风险表
           xlab = "Follow up time(M)",  #x轴标签
           break.x.by = 10,          #x轴刻度间距
           legend = c(0.8,0.75),     #图例位置
           legend.title = "",        #图例标题
) 



cor_data_df[which(cor_data_df$correlation>0),'Change'] <- 'Positive'
cor_data_df[which(cor_data_df$correlation<0),'Change'] <- 'Negative'





ggplot(data = cor_data_df, aes(x = correlation, y = -log10(fdr), color = Change)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c('red', 'blue'), limits = c('Positive', 'Negative')) +  #自定义点的颜色
  labs(x = 'correlation', y = '-log10 FDR', title = 'methlylation correlation', color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) 


cor_data_sig2 <- subset(cor_data_df, cor_data_df$fdr<0.05 & abs(cor_data_df$correlation) > 0.25)
table(cor_data_sig2$Change)
