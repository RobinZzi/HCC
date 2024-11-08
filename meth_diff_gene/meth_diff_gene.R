rm(list = ls())
library(data.table)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library("BioCor")
library(clusterProfiler)
library(enrichplot)
library(aPEAR)
library(msigdbr)
setwd("~/projects/hcc/analysis/meth_diff_gene")
save.image("meth_diff_gene.Rdata")


hcc2_uplist <- fread("pt1_pt2_up.txt")
hcc2_downlist <- fread("pt1_pt2_down.txt")


TCGA_pos_re <- fread("cor_data_df_pos_top.txt")
TCGA_neg_re <- fread("cor_data_df_neg_top.txt")

TCGA_diff <- fread("tcga_meth-type-based_wilcox-test_diff-gene_overlapwithhcc_result.txt")

TCGA_diff <- subset(TCGA_diff,subset=diff_fdr <0.05)
TCGA_diff_hypo_up <- subset(TCGA_diff,subset = log2FC > 0)
TCGA_diff_hypo_down<- subset(TCGA_diff,subset = log2FC < 0)

TCGA_diff_hypo_up <- subset(TCGA_diff,subset = log2FC > 0 & hypo_more_fdr <0.05)
TCGA_diff_hypo_down <- subset(TCGA_diff,subset = log2FC < 0& hypo_less_fdr <0.05)



intersect(intersect(TCGA_pos_re$symbol,hcc2_downlist$x),TCGA_diff_hypo_down$gene)
intersect(intersect(TCGA_neg_re$symbol,hcc2_uplist$x),TCGA_diff_hypo_up$gene)



up_venn <- list(TCGA_neg_re$symbol,hcc2_uplist$x,TCGA_diff_hypo_up$gene)
down_venn <- list(TCGA_pos_re$symbol,hcc2_downlist$x,TCGA_diff_hypo_down$gene)

venn.diagram(up_venn, filename = 'up.png', imagetype = 'png', ,category.names = c("tcga_neg" , "hypo_up" , "tcga_up"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')


venn.diagram(down_venn, filename = 'down.png', imagetype = 'png', ,category.names = c("tcga_pos" , "hypo_down" , "tcga_down"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')





hypo_up_go <- enrichGO(gene  = hcc1_uplist$x,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
hypo_up_go <- as.data.frame(hypo_up_go@result)
hypo_up_go [,"LogFDR"] <- -log10(hypo_up_go$p.adjust)
hypo_up_go [,"logp"] <- -log10(hypo_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
hypo_up_enrichment_fold <- apply(hypo_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
hypo_up_go$EF <- hypo_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = hypo_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Up-Regulated Pathways in de-methylated samples of HCC2")






hypo_down_go <- enrichGO(gene  = hcc1_downlist$x,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
hypo_down_go <- as.data.frame(hypo_down_go@result)
hypo_down_go [,"LogFDR"] <- -log10(hypo_down_go$p.adjust)
hypo_down_go [,"logp"] <- -log10(hypo_down_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
hyper_down_enrichment_fold <- apply(hypo_down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
hypo_down_go$EF <- hyper_down_enrichment_fold

#Draw Bubble Plot
ggplot(data = hypo_down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Down-Regulated Pathways in de-methylated samples of HCC2")



ggplot(data = hypo_down_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Down-Regulated Pathways in de-methylated samples")



methy_delta <- fread("methy_delta_screen.txt")


methy_delta_up <- subset(methy_delta,subset=DNA_meth_dt<0)
methy_delta_down <- subset(methy_delta,subset=DNA_meth_dt>0)

methy_delta_up_go <- enrichGO(gene  = methy_delta_up$gene,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
methy_delta_up_go <- as.data.frame(methy_delta_up_go@result)
methy_delta_up_go [,"LogFDR"] <- -log10(methy_delta_up_go$p.adjust)
methy_delta_up_go [,"logp"] <- -log10(methy_delta_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
methy_delta_up_enrichment_fold <- apply(methy_delta_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
methy_delta_up_go$EF <- methy_delta_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = methy_delta_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="de-methylated genes Pathways in samples of HCC2")




methy_delta_down_go <- enrichGO(gene  = methy_delta_down$gene,
                              OrgDb      = org.Hs.eg.db,
                              keyType    = 'SYMBOL',
                              ont        = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
methy_delta_down_go <- as.data.frame(methy_delta_down_go@result)
methy_delta_down_go [,"LogFDR"] <- -log10(methy_delta_down_go$p.adjust)
methy_delta_down_go [,"logp"] <- -log10(methy_delta_down_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
methy_delta_down_enrichment_fold <- apply(methy_delta_down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
methy_delta_down_go$EF <- methy_delta_down_enrichment_fold

#Draw Bubble Plot
ggplot(data = methy_delta_down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="ad-methylated genes Pathways in samples of HCC2")




trio_up <- intersect(methy_delta_up$gene,hcc1_uplist$x)

tcga_up <- intersect(TCGA_neg_re$symbol,TCGA_diff_hypo_up$gene)

tcga_down <- intersect(TCGA_pos_re$symbol,TCGA_diff_hypo_down$gene)

common_up <- intersect(trio_up,tcga_up)














tcga_up_go <- enrichGO(gene  = tcga_up,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
tcga_up_go <- as.data.frame(tcga_up_go@result)
tcga_up_go [,"LogFDR"] <- -log10(tcga_up_go$p.adjust)
tcga_up_go [,"logp"] <- -log10(tcga_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_up_enrichment_fold <- apply(tcga_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_up_go$EF <- tcga_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Up-Regulated Pathways in de-methylated samples of TCGA")



tcga_down_go <- enrichGO(gene  = tcga_down,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
tcga_down_go <- as.data.frame(tcga_down_go@result)
tcga_down_go [,"LogFDR"] <- -log10(tcga_down_go$p.adjust)
tcga_down_go [,"logp"] <- -log10(tcga_down_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_down_enrichment_fold <- apply(tcga_down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_down_go$EF <- tcga_down_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Down-Regulated Pathways in de-methylated samples of TCGA")



tcga_up_go_sig <- subset(tcga_up_go,pvalue < 0.05)












tcga_hypo_up_go <- enrichGO(gene  = TCGA_diff_hypo_up$gene,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
tcga_hypo_up_go <- as.data.frame(tcga_hypo_up_go@result)
tcga_hypo_up_go [,"LogFDR"] <- -log10(tcga_hypo_up_go$p.adjust)
tcga_hypo_up_go [,"logp"] <- -log10(tcga_hypo_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_up_enrichment_fold <- apply(tcga_hypo_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_hypo_up_go$EF <- tcga_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_hypo_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Up-Regulated Pathways in de-methylated samples of TCGA")


tcga_hypo_down_go <- enrichGO(gene  = TCGA_diff_hypo_down$gene,
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
tcga_hypo_down_go <- as.data.frame(tcga_hypo_down_go@result)
tcga_hypo_down_go [,"LogFDR"] <- -log10(tcga_hypo_down_go$p.adjust)
tcga_hypo_down_go [,"logp"] <- -log10(tcga_hypo_down_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_up_enrichment_fold <- apply(tcga_hypo_down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_hypo_down_go$EF <- tcga_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_hypo_down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Down-Regulated Pathways in de-methylated samples of TCGA")






tcga_neg_go <- enrichGO(gene  = TCGA_neg_re$symbol,
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
tcga_neg_go <- as.data.frame(tcga_neg_go@result)
tcga_neg_go [,"LogFDR"] <- -log10(tcga_neg_go$p.adjust)
tcga_neg_go [,"logp"] <- -log10(tcga_neg_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_neg_up_enrichment_fold <- apply(tcga_neg_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_neg_go$EF <- tcga_neg_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_neg_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Neg-Cor genes Pathways in samples of TCGA")




tcga_pos_go <- enrichGO(gene  = TCGA_pos_re$symbol,
                        OrgDb      = org.Hs.eg.db,
                        keyType    = 'SYMBOL',
                        ont        = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
tcga_pos_go <- as.data.frame(tcga_pos_go@result)
tcga_pos_go [,"LogFDR"] <- -log10(tcga_pos_go$p.adjust)
tcga_pos_go [,"logp"] <- -log10(tcga_pos_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_pos_enrichment_fold <- apply(tcga_pos_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_pos_go$EF <- tcga_pos_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_pos_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ",
       title="Pos-Cor genes Pathways in samples of TCGA")





summary_list <- list(methy_delta_up$gene,hcc1_uplist$x,TCGA_neg_re$symbol,TCGA_diff_hypo_up$gene)
venn.diagram(summary_list, filename = 'sum_up2.png', imagetype = 'png',category.names = c("hcc2_demethy" , "hcc2_rna_up" ,"tcga_methy_neg", "tcga_rna_up"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cat.cex = 0.8, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cex = 1.5, fontfamily = 'serif')




summary_down_list <- list(methy_delta_down$gene,hcc1_downlist$x,TCGA_pos_re$symbol,TCGA_diff_hypo_down$gene)
venn.diagram(summary_down_list, filename = 'sum_down2.png', imagetype = 'png', ,category.names = c("hcc2_admethy" , "hcc2_rna_down" ,"tcga_methy_pos", "tcga_rna_down"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cat.cex = 0.8, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cex = 1.5, fontfamily = 'serif')





methy_delta_up_go_sig <- subset(methy_delta_up_go,subset = p.adjust<0.05)
hypo_up_go_sig <- subset(hypo_up_go,subset = p.adjust<0.05)
tcga_neg_go_sig <- subset(tcga_neg_go,subset = p.adjust<0.05)
tcga_hypo_up_go_sig <-  subset(tcga_hypo_up_go,subset = p.adjust<0.05)

summary_up_pathway <- Reduce(rbind,list(methy_delta_up_go_sig,
                                            hypo_up_go_sig,
                                            tcga_neg_go_sig,
                                            tcga_hypo_up_go_sig))


methy_delta_up_go_sig <- subset(methy_delta_up_go,subset = p.adjust<0.05)
hypo_up_go_sig <- subset(hypo_up_go,subset = p.adjust<0.05)
tcga_neg_go_sig <- subset(tcga_neg_go,subset = p.adjust<0.05)
tcga_hypo_up_go_sig <-  subset(tcga_hypo_up_go,subset = p.adjust<0.05)

methy_delta_up_go_sig <- subset(methy_delta_up_go,subset = pvalue<0.05)
hypo_up_go_sig <- subset(hypo_up_go,subset = pvalue<0.05)
tcga_neg_go_sig <- subset(tcga_neg_go,subset = pvalue<0.05)
tcga_hypo_up_go_sig <-  subset(tcga_hypo_up_go,subset = pvalue<0.05)

summary_up_pathway <- Reduce(intersect,list(methy_delta_up_go_sig$Description,
                                            hypo_up_go_sig$Description,
                                            tcga_neg_go_sig$Description,
                                            tcga_hypo_up_go_sig$Description))


methy_delta_down_go_sig <- subset(methy_delta_down_go,subset = p.adjust<0.05)
hypo_down_go_sig <- subset(hypo_down_go,subset = p.adjust<0.05)
tcga_pos_go_sig <- subset(tcga_pos_go,subset = p.adjust<0.05)
tcga_hypo_down_go_sig <-  subset(tcga_hypo_down_go,subset = p.adjust<0.05)

summary_down_pathway <- Reduce(rbind,list(methy_delta_up_go_sig,
                                            hypo_up_go_sig,
                                            tcga_neg_go_sig,
                                            tcga_hypo_up_go_sig))









enrichmentNetwork(tcga_neg_go[1:200,],
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)

enrichmentNetwork(tcga_pos_go[1:200,],
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)

enrichmentNetwork(tcga_up_go[1:200,],
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)

enrichmentNetwork(tcga_down_go[1:200,],
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)

enrichmentNetwork(methy_delta_up_go_sig,
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)

enrichmentNetwork(methy_delta_down_go_sig,
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)


enrichmentNetwork(hypo_up_go_sig,
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)

enrichmentNetwork(hypo_down_go_sig,
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)





enrichmentNetwork(rbind(tcga_pos_go[1:50,],tcga_down_go[1:50,]),
                  colorBy='pvalue',
                  nodeSize='Count',
                  drawEllipses=TRUE)

summary_up_pathway <- summary_up_pathway[order(summary_up_pathway$p.adjust),]
enrichmentNetwork(summary_up_pathway[1:300,],
                  colorBy='p.adjust',
                  nodeSize='Count',
                  drawEllipses=TRUE,
                  allow.cartesian=TRUE)

summary_down_pathway <- summary_down_pathway[order(summary_down_pathway$p.adjust),]
enrichmentNetwork(summary_down_pathway[1:300,],
                  colorBy='p.adjust',
                  nodeSize='Count',
                  drawEllipses=TRUE,
                  allow.cartesian=TRUE)
