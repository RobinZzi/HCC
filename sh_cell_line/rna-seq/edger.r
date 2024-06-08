setwd("~/projects/hcc/analysis/rna-seq")
library(edgeR)
library(org.Hs.eg.db)
library(pheatmap)
library(stringr)
library(data.table)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(yulab.utils)

install.packages("Cairo")
library(Cairo)
options(bitmapType='cairo')
plot(cars)


rm(list=ls())
save.image("diff.Rdata")
mat = read.delim("combined-chrM.counts",header=T,skip=1)

counts = mat[,-c(1:6)]
counts = counts[,-c(5,12)]
rownames(counts) = mat$Geneid


group = c("SNU449-C","SNU449-C","SNU449-E","HepG2-E","HepG2-E","HepG2-C","HepG2-C","HepG2-C","SNU449-E","SNU449-E")
design = model.matrix(~0+group)


 y = DGEList(counts=counts)
 y<-calcNormFactors(y)
 y<-estimateCommonDisp(y, rowsum.filter=5)
 y<-estimateGLMTagwiseDisp(y,design) #
 pdf("MDSplot.pdf")
 plotMDS(y,label=group,xlim=c(-8,8))
 dev.off()
 
 
 fit_tag<-glmFit(y,design)
 lrt.tagwise<-glmLRT(fit_tag,contrast=c(1,-1,0,0))
 pvals_tag <- lrt.tagwise$table$PValue
 FDR_tag<- p.adjust(pvals_tag, method="BH")
 out1 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
 out1 = out1[order(out1$FDR_tag),]
 table(out1$FDR_tag<0.05)
 sig1 = out1[which(out1$FDR_tag<0.05),]
 setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq/diff_analysis")
 write.csv(sig1,"HepG2.diff.csv")

 lrt.tagwise<-glmLRT(fit_tag,contrast=c(0,0,1,-1))
 pvals_tag <- lrt.tagwise$table$PValue
 FDR_tag<- p.adjust(pvals_tag, method="BH")
 out2 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
 out2 = out2[order(out2$FDR_tag),]
 table(out2$FDR_tag<0.05)
 sig2 = out2[which(out2$FDR_tag<0.05),]
 setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq/diff_analysis")
 write.csv(sig2,"SNU449.diff.csv")

 
 lrt.tagwise<-glmLRT(fit_tag,contrast=c(1,0,-1,0))
 pvals_tag <- lrt.tagwise$table$PValue
 FDR_tag<- p.adjust(pvals_tag, method="BH")
 out3 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
 out3 = out3[order(out3$FDR_tag),]
 table(out3$FDR_tag<0.05)
 sig3 = out3[which(out3$FDR_tag<0.05),]
 setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq/diff_analysis")
 write.csv(sig3,"SNU449-vs-hip.diff.csv")

 
setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq/diff_analysis")
Hep_diff <- fread("HepG2.diff.csv")
Hep_up <- subset(Hep_diff,subset = logFC>0)
Hep_down <- subset(Hep_diff,subset = logFC<0)


Hep_up_go <- enrichGO(gene  = Hep_up$GeneID,
                               OrgDb      = org.Hs.eg.db,
                               keyType    = 'SYMBOL',
                               ont        = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05)
Hep_up_go <- as.data.frame(Hep_up_go@result)
Hep_up_go [,"logp"] <- -log10(Hep_up_go$pvalue)
ggplot(data = Hep_up_go[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Up-Regulated Pathways in sh-Gadd45a HepG2")



Hep_down_go <- enrichGO(gene  = Hep_down$GeneID,
                      OrgDb      = org.Hs.eg.db,
                      keyType    = 'SYMBOL',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
Hep_down_go <- as.data.frame(Hep_down_go@result)
Hep_down_go [,"logp"] <- -log10(Hep_down_go$pvalue)
ggplot(data = Hep_down_go[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Down-Regulated Pathways in sh-Gadd45a HepG2")

Hep_diff <- out1 %>%
  mutate(expression = case_when(logFC >= 1 & PValue < 0.05 ~ "Up-regulated", # 上调
                                logFC <= -1 & PValue < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变



ggplot(Hep_diff, aes(logFC, -log10(PValue))) +
  geom_point(size = 0.4, aes(color = expression)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" p-value")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="Hepg2 diff gene")


SNU_diff <- fread("SNU449.diff.csv")
SNU_up <- subset(SNU_diff,subset = logFC>0)
SNU_down <- subset(SNU_diff,subset = logFC<0)


SNU_up_go <- enrichGO(gene  = SNU_up$GeneID,
                      OrgDb      = org.Hs.eg.db,
                      keyType    = 'SYMBOL',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
SNU_up_go <- as.data.frame(SNU_up_go@result)
SNU_up_go [,"logp"] <- -log10(SNU_up_go$pvalue)
ggplot(data = SNU_up_go[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Up-Regulated Pathways in sh-Gadd45a SNU449")



SNU_down_go <- enrichGO(gene  = SNU_down$GeneID,
                        OrgDb      = org.Hs.eg.db,
                        keyType    = 'SYMBOL',
                        ont        = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
SNU_down_go <- as.data.frame(SNU_down_go@result)
SNU_down_go [,"logp"] <- -log10(SNU_down_go$pvalue)
ggplot(data = SNU_down_go[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Down Pathways in sh-Gadd45a SNU449")


SNU_diff <- out2 %>%
  mutate(expression = case_when(logFC >= 1 & PValue < 0.05 ~ "Up-regulated", # 上调
                                logFC <= -1 & PValue < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变



ggplot(SNU_diff, aes(logFC, -log10(PValue))) +
  geom_point(size = 0.4, aes(color = expression)) +
  xlab(expression("log"[2]*" fold change")) + 
  ylab(expression("-log"[10]*" p-value")) + 
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ 
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="SNU449 diff gene")






SH_diff <- fread("SNU449-vs-hEp.diff.csv")
SH_up <- subset(SH_diff,subset = logFC>0)
SH_down <- subset(SH_diff,subset = logFC<0)


SH_up_go <- enrichGO(gene  = SH_up$GeneID,
                      OrgDb      = org.Hs.eg.db,
                      keyType    = 'SYMBOL',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
SH_up_go <- as.data.frame(SH_up_go@result)
SH_up_go [,"logp"] <- -log10(SH_up_go$pvalue)

SH_up_go_ch <- mutate(SH_up_go,Description = en2ch(Description))
ggplot(data = SH_up_go[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Up-Regulated Pathways in sh-Gadd45a SNU449")



SH_down_go <- enrichGO(gene  = SH_down$GeneID,
                        OrgDb      = org.Hs.eg.db,
                        keyType    = 'SYMBOL',
                        ont        = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
SH_down_go <- as.data.frame(SH_down_go@result)
SH_down_go [,"logp"] <- -log10(SH_down_go$pvalue)
ggplot(data = SH_down_go[1:20,])+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Down Pathways in sh-Gadd45a SNU449")


SH_diff_all <- out3 %>%
  mutate(expression = case_when(logFC >= 1 & PValue < 0.05 ~ "Up-regulated", # 上调
                                logFC <= -1 & PValue < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变



ggplot(SH_diff_all, aes(logFC, -log10(PValue))) +
  geom_point(size = 0.4, aes(color = expression)) +
  xlab(expression("log"[2]*" fold change")) + 
  ylab(expression("-log"[10]*" p-value")) + 
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ 
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="SNU449 diff gene")




hpc_up <- fread("hpc_up_markers.txt")
hpc_down <- fread("hpc_down_markers.txt")



up_intersect_hep <- intersect(hpc_up$V2,Hep_down$GeneID)
down_intersect_hep<- intersect(hpc_down$V2,Hep_up$GeneID)

up_intersect_snu <- intersect(hpc_up$V2,SNU_down$GeneID)
down_intersect_snu <- intersect(hpc_down$V2,SNU_up$GeneID)


require(clusterProfiler)
library(devtools)

install_github("YuLab-SMU/yulab.utils")
1


