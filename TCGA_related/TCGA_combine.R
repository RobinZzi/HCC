

setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/TCGA-DNA_methylation/LIHC")
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
  library(BiocManager)
  library(ggpubr)
  library(tidyr)
  library(survminer) 
  library(ggcorrplot)
  library(survival)
  library(edgeR)
  library(R.utils)
  library(stringr)
  library(tinyarray)
  library(devtools)
}
rm(list=ls())

load("TCGA-LIHC.Rdata")
save.image("TCGA-LIHC.Rdata")

meth <- fread("TCGA-LIHC.methylation450.tsv",header = T)
pheno <- fread("TCGA-LIHC.phenotype.tsv",header = T)
anno <- fread("illuminaMethyl450_hg38_GDC")
fpkm <- fread("TCGA-LIHC.htseq_fpkm.tsv")
surv <- fread("TCGA-LIHC.survival.tsv")
count <- fread("TCGA-LIHC.htseq_counts.tsv.gz")
cnv <- fread("TCGA-LIHC.cnv.tsv",header = T)
mut <- fread("TCGA-LIHC.mutect2_snv.tsv",header = T)

#meth_filt#
meth_filt <-subset(meth,meth$`Composite Element REF` %in% anno$`#id`)
rownames(meth_filt) <- anno$`#id`
meth_filt2 <- na.omit(meth_filt)
meth_filt2 <- meth_filt2[,-1]


rank <- fread("rank.txt")
write.csv(rank,"rank.csv")



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
               y = XRN2,
               x = mean_meth,
               centrality.para = "mean"
)


ggscatterstats(data = comb, 
               y = GADD45A,
               x = mean_meth,
               centrality.para = "mean"
)





VlnPlot(hpc,features = "UHRF1",group.by = "sample")


exp_mtx <- fread("exp_meanmeth_mtx.txt")
exp_mtx <- exp_mtx[,-2]
exp_mtx_t <-as.data.frame(t(exp_mtx))
sampleid_mtx <- exp_mtx_t[1,]
sampleid_mtx  <- unlist(sampleid_mtx)
colnames(exp_mtx_t) <- sampleid_mtx
exp_mtx_t <- exp_mtx_t[-1,]


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
               y = SETDB1, 
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

pheatmap(sort_TAGA_data2,
         show_rownames = F, show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 100,
         treeheight_col = 20,
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




#id-transform#

Ensembl_ID <- count_foridtransfer$Ensembl_ID
Ensembl_ID <- str_sub(Ensembl_ID,1,15)
count_foridtransfer$Ensembl_ID2 <- Ensembl_ID
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
gene_symbol <- fread("gene_symbol.txt")
gene_symbol <- gene_symbol[,-1]
count_filt <- subset(count_foridtransfer, count_foridtransfer$Ensembl_ID2 %in% gene_symbol$ENSEMBL)
gene_symbol <- gene_symbol[match(count_filt$Ensembl_ID2,gene_symbol$ENSEMBL),]
count_filt$SYMBOL <- gene_symbol$SYMBOL
exp <- count_filt
exp2<- avereps(exp,ID=exp$SYMBOL)
exp2 <- as.data.frame(exp2)
exp2 <- exp2[!is.na(exp2$SYMBOL),]#去除没有匹配的
exp2 <- exp2[,-1]
row.names(exp2) <- exp2$SYMBOL
exp3 <- select(exp2,!SYMBOL)
exp3 <- select(exp3,count_order_colnames)
exp3 <- as.data.frame(lapply(exp3,as.numeric))
row.names(exp3) <- exp2$SYMBOL

  dge <- DGEList(counts=exp3,group=Group)
  
  dge$samples$lib.size <- colSums(dge$counts)
  dge <- calcNormFactors(dge) 
  
  design <- model.matrix(~0+Group)
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(Group)
  
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  fit2 <- glmLRT(fit, contrast=c(-1,1)) 
  
  DEG=topTags(fit2, n=nrow(exp))
  DEG=as.data.frame(DEG)
  
  logFC_t=0
  P.Value_t = 0.05
  k1 = (DEG$FDR < P.Value_t)&(DEG$logFC < -logFC_t)
  k2 = (DEG$FDR < P.Value_t)&(DEG$logFC > logFC_t)
  change = ifelse(k1,"down",ifelse(k2,"up","stable"))
  DEG$change <- change
  
  
  DEG_edgeR <- DEG
  DEG_edgeR[,'gene'] <- row.names(DEG_edgeR)
  
  table(DEG_edgeR$change)


write.table(DEG_edgeR,"DEG_edgeR_hypoall_vs_hyper.txt")

DEG_edgeR_hypoup <- subset(DEG_edgeR, subset = DEG_edgeR$change == 'up')
DEG_edgeR_hypodown <- subset(DEG_edgeR, subset = DEG_edgeR$change == 'down')


write.table(hypo_up_bp,"DEG_hypoup_bp_hypoall_vs_hyper.txt")
write.table(hypo_down_bp,"DEG_hypodown_bp_hypoall_vs_hyper.txt")


hcc3_demeth_down <- fread("hcc3_demeth_down.txt")



hcc_markers_hypoup <- subset(hcc_markers, subset = hcc_markers$avg_log2FC > 0.5 & hcc_markers$p_val_adj < 0.05)
hcc_markers_hypodown <- subset(hcc_markers, subset = hcc_markers$avg_log2FC < -0.5 & hcc_markers$p_val_adj < 0.05)




hcc_markers_commonup <- subset(hcc_markers_hypoup, subset= hcc_markers_hypoup$gene %in% DEG_edgeR_hypodown$gene)
hcc_markers_commondown <- subset(hcc_markers_hypoup, subset= hcc_markers_hypoup$gene %in% DEG_edgeR_hypodown$gene)







count_order_tr <- as.data.frame(t(count_order_data))

count_order_tr[,"sample_type"] <- as.factor(c(rep("hypo",261),rep("hyper",150)))

count_order_tr <- select(count_order_tr,"sample_type",everything())
ggplot(cpm_order_tr) +
  aes(x = sample_type, y = GADD45A) +
  geom_boxplot(fill = "#0c4c8a") +
  theme_minimal()
hist(subset(cpm_order_tr, sample_type == "hypo")$GADD45A,
     main = "hypo expression",
     xlab = "expression"
)
hist(subset(cpm_order_tr, sample_type == "hyper")$GADD45A,
     main = "hyper expression",
     xlab = "expression"
)
shapiro.test(subset(cpm_order_tr, sample_type == "hypo")$GADD45A)
shapiro.test(subset(cpm_order_tr, sample_type == "hyper")$GADD45A)
wilcox.test(cpm_order_tr$GADD45A ~ cpm_order_tr$sample_type)

cpm_order_data <- edgeR::cpm(count_order_data)
cpm_order_tr <- as.data.frame(t(cpm_order_data))
cpm_order_tr[,"sample_type"] <- as.factor(c(rep("hypo",261),rep("hyper",150)))
cpm_order_tr <- select(cpm_order_tr,"sample_type",everything())

lihc_meth  <- numeric(6)
lihc_methdiff_base <- as.data.frame(lihc_meth)
cpm_order_tr_hyper <- subset(cpm_order_tr , subset = cpm_order_tr$sample_type == "hyper")
cpm_order_tr_hypo <- subset(cpm_order_tr , subset = cpm_order_tr$sample_type == "hypo")
for(i in 2: length(cpm_order_tr)){
  test_diff <- wilcox.test(cpm_order_tr[,i] ~ cpm_order_tr$sample_type)
  test_diff_less <- wilcox.test(cpm_order_tr[,i] ~ cpm_order_tr$sample_type,alternative = "less")
  test_diff_more <- wilcox.test(cpm_order_tr[,i] ~ cpm_order_tr$sample_type,alternative = "greater")
  hypomean <- mean(log2(cpm_order_tr_hypo[,names(cpm_order_tr)[i]]+1))
  hypermean <- mean(log2(cpm_order_tr_hyper[,names(cpm_order_tr)[i]]+1))
  logFC <- hypomean-hypermean
  lihc_meth[1] <- test_diff$p.value
  lihc_meth[2] <- test_diff_less$p.value
  lihc_meth[3] <- test_diff_more$p.value
  lihc_meth[4] <- hypomean
  lihc_meth[5] <- hypermean
  lihc_meth[6] <- logFC
  lihc_methdiff_base <- cbind(lihc_methdiff_base,lihc_meth)
  names(lihc_methdiff_base)[i] <- names(cpm_order_tr)[i]
}

lihc_methdiff_result <- as.data.frame(t(lihc_methdiff_base))
colnames(lihc_methdiff_result) <- c("diff","hypo_more","hypo_less","hypomean","hypermean","log2FC")
lihc_methdiff_result <- lihc_methdiff_result[-1,]
lihc_methdiff_result$gene <- row.names(lihc_methdiff_result)
lihc_methdiff_result$diff_padj <- p.adjust(as.vector(lihc_methdiff_result$diff), "fdr", n = length(lihc_methdiff_result$diff))
lihc_methdiff_result$hypo_more_padj <- p.adjust(as.vector(lihc_methdiff_result$hypo_more), "fdr", n = length(lihc_methdiff_result$hypo_more))
lihc_methdiff_result$hypo_less_padj <- p.adjust(as.vector(lihc_methdiff_result$hypo_less), "fdr", n = length(lihc_methdiff_result$hypo_less))

lihc_methdiff_result$change <- ifelse(lihc_methdiff_result$hypo_more_padj < 0.05,"up",ifelse(lihc_methdiff_result$hypo_less_padj < 0.05,"down","stable"))
lihc_methdiff_result$change2 <- ifelse(lihc_methdiff_result$diff_padj < 0.05,"change","stable")





lihc_methdiff_result_order <- lihc_methdiff_result[order(lihc_methdiff_result$diff_padj),]
lihc_methdiff_result_order$diff_order <- c(1:nrow(lihc_methdiff_result_order))


colnames(lihc_methdiff_result_order) <- c("diff_p.value","hypo_more_p.value","hypo_less_p.value","hypo_mean_logCPM","hyper_mean_logCPM","log2FC","gene","diff_fdr","hypo_more_fdr","hypo_less_fdr","change_type","state","diff_fdr_rank")


lihc_methdiff_result_order <- select(lihc_methdiff_result_order,c("gene","hypo_mean_logCPM","hyper_mean_logCPM","log2FC","diff_p.value","diff_fdr","hypo_more_p.value","hypo_more_fdr","hypo_less_p.value","hypo_less_fdr","state","change_type","diff_fdr_rank"))

write.table(lihc_methdiff_result_order,"tcga_meth-type-based_wilcox-test_diff-gene_result.txt")









single_cell_deg_result_up <- fread("hpc_up_markers.txt")
single_cell_deg_result_down <- fread("hpc_down_markers.txt")


lihc_methdiff_result_order_up <- subset(lihc_methdiff_result_order, subset = lihc_methdiff_result_order$change_type == 'up')
lihc_methdiff_result_order_down <- subset(lihc_methdiff_result_order, subset = lihc_methdiff_result_order$change_type == 'down')
common_up <- intersect(single_cell_deg_result_up$V2,lihc_methdiff_result_order_up$gene)
common_down <- intersect(single_cell_deg_result_down$V2,lihc_methdiff_result_order_down$gene)







hypo_up_bp <-enrichGO(gene       = common_up,
                      OrgDb      = org.Hs.eg.db,
                      keyType    = 'SYMBOL',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
hypo_up_bp <- as.data.frame(hypo_up_bp@result)
hypo_up_bp[,"logp"] <- -log10(hypo_up_bp$pvalue)
ggplot(data = hypo_up_bp[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  labs(x = "-Logp", 
       y= " ")






hypo_down_bp <-enrichGO(gene       = common_down,
                        OrgDb      = org.Hs.eg.db,
                        keyType    = 'SYMBOL',
                        ont        = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
hypo_down_bp <- as.data.frame(hypo_down_bp@result)
hypo_down_bp[,"logp"] <- -log10(hypo_down_bp$pvalue)
ggplot(data = hypo_down_bp[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  labs(x = "-Logp", 
       y= " ")




overlap_result <- rbind(lihc_methdiff_result_order_up,lihc_methdiff_result_order_down)
write.table(overlap_result,"tcga_meth-type-based_wilcox-test_diff-gene_overlapwithhcc_result.txt")









sub_mut_info_order <-  sub_mut_info[order(sub_mut_info$fisher_pvalue),]


filt_sub_mut_info <- subset(sub_mut_info_order, subset = sub_mut_info_order$fisher_pvalue < 0.05)
filt_sub_mut_info_hypo <- subset(filt_sub_mut_info,subset= filt_sub_mut_info$enrich_factor == "hypo")
filt_sub_mut_info_hyper <- subset(filt_sub_mut_info,subset= filt_sub_mut_info$enrich_factor == "hyper")
write.table(sub_mut_info_order,"TCGA-LIHC_methgroup_fisher_test_result.txt")
write.table(filt_sub_mut_info,"TCGA-LIHC_methgroup_fisher_test_result_filtbyPvalue.txt")
write.table(filt_sub_mut_info_hypo,"TCGA-LIHC_methgroup_fisher_test_result_filtbyPvalue_hypo.txt")
write.table(filt_sub_mut_info_hyper,"TCGA-LIHC_methgroup_fisher_test_result_filtbyPvalue_hyper.txt")



write.table(sort_TAGA_data2,"TCGA-LIHC_sort_meth_level.txt")












sort_TAGA_data = data.frame(test_merge[TCGA_order_row,TCGA_order_col])




baseMean = colMeans(sort_TAGA_data2)

pheatmap(baseMean,show_col_names = F,cluster_rows = F, cluster_cols = F)
colanno <- colnames(sort_TAGA_data2)
colanno <- as.data.frame(colanno)
colanno_plot$sample_name = factor(colanno_plot$sample_name, levels = colnames(sort_TAGA_data2))



colanno$mean <- baseMean
row.names(colanno) <- colanno$colanno
colnames(colanno) <- c("sample_name","mean")
colanno_plot <- colanno


pheatmap(sort_TAGA_data2,
         show_rownames = F, show_colnames = F,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 100,
         treeheight_col = 20,
         fontsize          = 12)


ggplot(colanno_plot,aes(x=sample_name,y=mean,fill=mean))+
  geom_point(size=0.7,shape = 21, stroke = 0)+
  scale_fill_steps(low = "white", high = "red")+theme_classic()


ggplot(colanno_plot,aes(x=sample_name,y=mean))+
  geom_point(size=0.7,shape = 16)+theme_classic()


dev.off()
