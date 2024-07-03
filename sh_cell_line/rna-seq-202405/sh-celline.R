rm(list = ls())
library(FactoMineR)
setwd("~/projects/hcc/analysis/liver_cancer_cellline")
count2 <- fread("combined-chrM.counts")
rpkm2 <- fread("combined-chrM.rpkm")

setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-202405")

count <- fread("combined-chrM.counts")
rpkm <- fread("combined-chrM.rpkm")

setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-202405/diff_analysis")
save.image("sh-celline.Rdata")








sum_rpkm <- as.data.frame(cbind(select(rpkm,1,7,8,9,10),rpkm2[,7]))
rownames(sum_rpkm) <- sum_rpkm[,1]
sum_rpkm <- sum_rpkm[,-1]
colnames(sum_rpkm) <- c("sh-1","sh-2","v-1","v-2","wt")



sum_count <- as.data.frame(cbind(select(count,1,7,8,9,10),count2[,7]))
rownames(sum_count) <- sum_count[,1]
sum_count <- sum_count[,-1]
colnames(sum_count) <- c("sh-1","sh-2","v-1","v-2","wt")


genelist <- c("SNHG6","MAT1A","MAT2A","FUS")

submtx <- sum_rpkm[genelist,]
submtx <- t(submtx)


group = c("sh","sh","v","v","wt")
design = model.matrix(~0+group)


y = DGEList(counts=sum_count)
y<-calcNormFactors(y)
y<-estimateCommonDisp(y, rowsum.filter=5)
y<-estimateGLMTagwiseDisp(y,design) #
pdf("MDSplot.pdf")
plotMDS(y,label=group,xlim=c(-8,8))
dev.off()


fit_tag<-glmFit(y,design)
lrt.tagwise<-glmLRT(fit_tag,contrast=c(1,-1,0))
pvals_tag <- lrt.tagwise$table$PValue
FDR_tag<- p.adjust(pvals_tag, method="BH")
out1 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
out1 = out1[order(out1$FDR_tag),]
table(out1$FDR_tag<0.05)
sig1 = out1[which(out1$FDR_tag<0.05),]
setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-202405/diff_analysis")
write.csv(sig1,"sh-v-diff.csv")

lrt.tagwise<-glmLRT(fit_tag,contrast=c(0,1,-1))
pvals_tag <- lrt.tagwise$table$PValue
FDR_tag<- p.adjust(pvals_tag, method="BH")
out2 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
out2 = out2[order(out2$FDR_tag),]
table(out2$FDR_tag<0.05)
sig2 = out2[which(out2$FDR_tag<0.05),]
setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-202405/diff_analysis")
write.csv(sig2,"v-wt-diff.csv")


lrt.tagwise<-glmLRT(fit_tag,contrast=c(1,0,-1))
pvals_tag <- lrt.tagwise$table$PValue
FDR_tag<- p.adjust(pvals_tag, method="BH")
out3 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
out3 = out3[order(out3$FDR_tag),]
table(out3$FDR_tag<0.05)
sig3 = out3[which(out3$FDR_tag<0.05),]
setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-202405/diff_analysis")
write.csv(sig3,"sh-wt.diff.csv")





setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-202405/diff_analysis")
sh_v_diff <- fread("sh-v-diff.csv")
sh_v_up <- subset(sh_v_diff,subset = logFC>0)
sh_v_down <- subset(sh_v_diff,subset = logFC<0)

sh_v_up_go <- enrichGO(gene  = sh_v_up$GeneID,
                      OrgDb      = org.Hs.eg.db,
                      keyType    = 'SYMBOL',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
sh_v_up_go <- as.data.frame(sh_v_up_go@result)
sh_v_up_go [,"logp"] <- -log10(sh_v_up_go$pvalue)
sh_v_up_go[,11:12] <- as.numeric(str_split_fixed(sh_v_up_go$GeneRatio,"/",2))
sh_v_up_go$GeneRatio <- sh_v_up_go[,11]/sh_v_up_go[,12]
ggplot(data = sh_v_up_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Up-Regulated Pathways in sh-SNHG6 Huh7")



sh_v_down_go <- enrichGO(gene  = sh_v_down$GeneID,
                         OrgDb      = org.Hs.eg.db,
                         keyType    = 'SYMBOL',
                         ont        = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
sh_v_down_go <- as.data.frame(sh_v_down_go@result)
sh_v_down_go [,"logp"] <- -log10(sh_v_down_go$pvalue)
sh_v_down_go[,11:12] <- as.numeric(str_split_fixed(sh_v_down_go$GeneRatio,"/",2))
sh_v_down_go$GeneRatio <- sh_v_down_go[,11]/sh_v_down_go[,12]
ggplot(data = sh_v_down_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="down-Regulated Pathways in sh-SNHG6 Huh7")

sh_v_diff_all <- out1 %>%
  mutate(expression = case_when(logFC >= 1 & FDR_tag < 0.05 ~ "Up-regulated", # 上调
                                logFC <= -1 & FDR_tag < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变



ggplot(sh_v_diff_all, aes(logFC, -log10(FDR_tag))) +
  geom_point(size = 0.4, aes(color = expression)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="sh-SNHG6 vs Vector diff gene")





v_wt_diff_all <- out2 %>%
  mutate(expression = case_when(logFC >= 1 & FDR_tag < 0.05 ~ "Up-regulated", # 上调
                                logFC <= -1 & FDR_tag < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变



ggplot(v_wt_diff_all, aes(logFC, -log10(FDR_tag))) +
  geom_point(size = 0.4, aes(color = expression)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Vector vs WT diff gene")


sh_wt_diff_all <- out3 %>%
  mutate(expression = case_when(logFC >= 1 & FDR_tag < 0.05 ~ "Up-regulated", # 上调
                                logFC <= -1 & FDR_tag < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变



ggplot(sh_wt_diff_all, aes(logFC, -log10(FDR_tag))) +
  geom_point(size = 0.4, aes(color = expression)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="sh-SNHG6 vs WT diff gene")


example <- sum_rpkm[c("ANO1","EHD2","UBE2E2"),]
example2 <- subset(out1,subset = GeneID %in% c("ANO1","EHD2","UBE2E2"))
example <- cbind(select(example2,1,2,6),example)


rpkm_t <- t(sum_rpkm)
rpkm_group <- data.frame(Sample = rownames(rpkm_t), Group = c("sh","sh","v","v","wt"))
rpkm_pca <- PCA(rpkm_t,scale.unit = TRUE,graph = FALSE)
rpkm_pca_sample <- data.frame(rpkm_pca$ind$coord[ ,1:2])
rpkm_pca_sample$Sample=row.names(rpkm_pca_sample)
rpkm_pca_eig1 <- round(rpkm_pca$eig[1,2], 2)
rpkm_pca_eig2 <- round(rpkm_pca$eig[2,2],2 )
rpkm_pca_sample <- merge(rpkm_pca_sample,rpkm_group,by="Sample")
head(rpkm_pca_sample)
ggplot(data = rpkm_pca_sample, aes(x = Dim.1, y = Dim.2))+
  geom_point(aes(color = Group), size = 5) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', rpkm_pca_eig1, '%'), y = paste('PCA2:', rpkm_pca_eig2, '%'), color = '') + 
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="SNHG6_pca")

out_list <- out1[genelist,]
out_list <- dplyr::select(out_list,1,2,6)
