setwd("~/projects/hcc/analysis/rna-seq")
mat = read.delim("combined-chrM.counts",header=T,skip=1)

counts = mat[,-c(1:6)]
counts = counts[,-c(2,10)]
rownames(counts) = mat$Geneid
group = sub("bam.*(shG45A.*).._combined.sorted.bam","\\1",colnames(counts))
design = model.matrix(~0+group)

library(edgeR)
library(org.Hs.eg.db)
library(pheatmap)
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
write.csv(sig1,"HepG2.diff.csv")

 lrt.tagwise<-glmLRT(fit_tag,contrast=c(0,0,1,-1))
 pvals_tag <- lrt.tagwise$table$PValue
 FDR_tag<- p.adjust(pvals_tag, method="BH")
 out2 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
 out2 = out2[order(out2$FDR_tag),]
table(out2$FDR_tag<0.05)
 sig2 = out2[which(out2$FDR_tag<0.05),]
write.csv(sig2,"SNU449.diff.csv")



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



phenoData()