rm(list=ls())

setwd("~/projects/hcc/analysis/meth_diff_gene/intersect")
load("intersect.Rdata")
save.image("intersect.Rdata")
methy_delta <- fread("methy_change.txt")


methy_delta_de <- subset(methy_delta,subset=dt<0)
methy_delta_ad <- subset(methy_delta,subset=dt>0)
cdgenes <- c("CD1A","CD1B","CD1C","CD1D","CD1E","FI16","AIM2","PYHIN1","MNDA")

methy_delta_cd <- subset(methy_delta,subset=gene %in% cdgenes)
cdedge <- subset(tcga_edgeR, subset=gene %in% cdgenes)
cdwilcox <- subset()
pt_up <- fread("pt1_pt2_up.txt")
pt_down <- fread("pt1_pt2_down.txt")


hcc2_up <- fread("hcc2_hypo_up_genes.txt")


tcga_edgeR <- fread("DEG_edgeR_hypo_vs_hyper.txt")
tcga_edgeR_up <- subset(tcga_edgeR,subset=change=="up")
tcga_edgeR_down <- subset(tcga_edgeR,subset=change=="down")

tcga_edgeR_fc <- fread("DEG_edgeR_hypo_vs_hyper_fc1.txt")
tcga_edgeR_fc_up <- subset(tcga_edgeR_fc,subset=change=="up")
tcga_edgeR_fc_down <- subset(tcga_edgeR_fc,subset=change=="down")

tcga_wilcox <- fread("tcga_meth-type-based_wilcox-test_diff-gene_result_hypo1vhyper.txt")
tcga_wilcox_up <- subset(tcga_wilcox,subset=change=="up")
tcga_wilcox_down <- subset(tcga_wilcox,subset=change=="down")




tcga_cor_pos <- fread("cor_data_df_pos_top.csv")
tcga_cor_pos <- tcga_cor_pos[2:1001,]
tcga_cor_neg <- fread("cor_data_df_neg_top.csv")

demethy_list <- list(hcc2_up$gene,tcga_cor_neg$symbol,tcga_wilcox_up$gene)

demethy_list <- list(pt_up$x,tcga_cor_neg$symbol,tcga_wilcox_up$gene)
demethy_sum <- Reduce(intersect,demethy_list)
venn.diagram(demethy_list, filename = 'up_de_0818.png', imagetype = 'png',
             category.names = c("hypo_up" ,"tcga_cor_neg","tcga_wilcox_up"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 0.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

admethy_list <- list(methy_delta_ad$gene,tcga_edgeR_down$gene,tcga_wilcox_down$gene,tcga_cor_pos$symbol)
admethy_sum <- Reduce(intersect,admethy_list)
venn.diagram(admethy_list, filename = 'down_ad.png', imagetype = 'png',
             category.names = c("hypo_down" , "hypo_admeth" , "tcga_edgeR_down","tcga_wilcox_down","tcga_cor_pos"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), cat.cex =0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), cex = 1, fontfamily = 'serif')
















hypo_up_go <- enrichGO(gene  = pt_up$x,
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









hypo_down_go <- enrichGO(gene  = pt_down$x,
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











methy_delta_up_go <- enrichGO(gene  = methy_delta_ad$gene,
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




methy_delta_down_go <- enrichGO(gene  = methy_delta_de$gene,
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
       title="de-methylated genes Pathways in samples of HCC2")




methy_delta_up_go <- enrichGO(gene  = methy_delta_ad$gene,
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
       title="ad-methylated genes Pathways in samples of HCC2")




tcga_edgeR_up_go <- enrichGO(gene  = tcga_edgeR_up$gene,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
tcga_edgeR_up_go <- as.data.frame(tcga_edgeR_up_go@result)
tcga_edgeR_up_go [,"LogFDR"] <- -log10(tcga_edgeR_up_go$p.adjust)
tcga_edgeR_up_go [,"logp"] <- -log10(tcga_edgeR_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_up_enrichment_fold <- apply(tcga_edgeR_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_edgeR_up_go$EF <- tcga_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_edgeR_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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



tcga_edgeR_down_go <- enrichGO(gene  = tcga_edgeR_down$gene,
                         OrgDb      = org.Hs.eg.db,
                         keyType    = 'SYMBOL',
                         ont        = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
tcga_edgeR_down_go <- as.data.frame(tcga_edgeR_down_go@result)
tcga_edgeR_down_go [,"LogFDR"] <- -log10(tcga_edgeR_down_go$p.adjust)
tcga_edgeR_down_go [,"logp"] <- -log10(tcga_edgeR_down_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_down_enrichment_fold <- apply(tcga_edgeR_down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_edgeR_down_go$EF <- tcga_down_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_edgeR_down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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







tcga_wilcox_up_go <- enrichGO(gene  = tcga_wilcox_up$gene,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
tcga_wilcox_up_go <- as.data.frame(tcga_wilcox_up_go@result)
tcga_wilcox_up_go [,"LogFDR"] <- -log10(tcga_wilcox_up_go$p.adjust)
tcga_wilcox_up_go [,"logp"] <- -log10(tcga_wilcox_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_up_enrichment_fold <- apply(tcga_wilcox_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_wilcox_up_go$EF <- tcga_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_wilcox_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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



tcga_wilcox_down_go <- enrichGO(gene  = tcga_wilcox_down$gene,
                         OrgDb      = org.Hs.eg.db,
                         keyType    = 'SYMBOL',
                         ont        = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
tcga_wilcox_down_go <- as.data.frame(tcga_wilcox_down_go@result)
tcga_wilcox_down_go [,"LogFDR"] <- -log10(tcga_wilcox_down_go$p.adjust)
tcga_wilcox_down_go [,"logp"] <- -log10(tcga_wilcox_down_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_down_enrichment_fold <- apply(tcga_wilcox_down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_wilcox_down_go$EF <- tcga_down_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_wilcox_down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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










tcga_neg_go <- enrichGO(gene  = tcga_cor_neg$symbol,
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




tcga_pos_go <- enrichGO(gene  = tcga_cor_pos$symbol,
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

tcga_edgeR_fc_up_go <- enrichGO(gene  = tcga_edgeR_fc_up$gene,
                                OrgDb      = org.Hs.eg.db,
                                keyType    = 'SYMBOL',
                                ont        = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
tcga_edgeR_fc_up_go <- as.data.frame(tcga_edgeR_fc_up_go@result)
tcga_edgeR_fc_up_go [,"LogFDR"] <- -log10(tcga_edgeR_fc_up_go$p.adjust)
tcga_edgeR_fc_up_go [,"logp"] <- -log10(tcga_edgeR_fc_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_up_enrichment_fold <- apply(tcga_edgeR_fc_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_edgeR_fc_up_go$EF <- tcga_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_edgeR_fc_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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
       title="Up-Regulated Pathways in de-methylated samples of TCGA based on edgeR (|log2FC|>1,fdr<0.05|)")



tcga_edgeR_fc_down_go <- enrichGO(gene  = tcga_edgeR_fc_down$gene,
                                  OrgDb      = org.Hs.eg.db,
                                  keyType    = 'SYMBOL',
                                  ont        = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
tcga_edgeR_fc_down_go <- as.data.frame(tcga_edgeR_fc_down_go@result)
tcga_edgeR_fc_down_go [,"LogFDR"] <- -log10(tcga_edgeR_fc_down_go$p.adjust)
tcga_edgeR_fc_down_go [,"logp"] <- -log10(tcga_edgeR_fc_down_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tcga_down_enrichment_fold <- apply(tcga_edgeR_fc_down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tcga_edgeR_fc_down_go$EF <- tcga_down_enrichment_fold

#Draw Bubble Plot
ggplot(data = tcga_edgeR_fc_down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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
       title="Down-Regulated Pathways in de-methylated samples of TCGA(|log2FC|>1,fdr<0.05|)")


cr_v_diff_all <- out1 %>%
  mutate(expression = case_when(logFC <= -1 & FDR_tag < 0.05 ~ "Up-regulated", # 上调
                                logFC >= 1 & FDR_tag < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变


write.table(tcga_edgeR_up_go,"tcga_edgeR_up_go.csv")
write.table(tcga_edgeR_up_go,"tcga_edgeR_up_go.csv")

ggplot(tcga_edgeR, aes(logFC, -log10(FDR))) +
  geom_point(size = 0.4, aes(color = matabolic_gene)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-5, 5)) +
  scale_color_manual(values = c( "red","grey"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="tcga_edgeR diff gene")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey")

ggplot(tcga_wilcox, aes(log2FC, -log10(diff_padj))) +
  geom_point(size = 0.4, aes(color = change)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-5, 5)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="tcga_wilcox diff gene")






methy_delta_up_go_sig <- subset(methy_delta_up_go,subset = p.adjust<0.05)
hypo_up_go_sig <- subset(hypo_up_go,subset = p.adjust<0.05)
tcga_neg_go_sig <- subset(tcga_neg_go,subset = p.adjust<0.05)
tcga_edgeR_up_go_sig <-  subset(tcga_edgeR_up_go,subset = p.adjust<0.05)
tcga_wilcox_up_go_sig <-  subset(tcga_wilcox_up_go,subset = p.adjust<0.05)

summary_up_pathway <- list(hypo_up_go_sig$Description,
                           methy_delta_up_go_sig$Description,
                           tcga_edgeR_up_go_sig$Description,
                           tcga_wilcox_up_go_sig$Description,
                           tcga_neg_go_sig$Description)

intersection_up_pathway <- Reduce(intersect,summary_up_pathway)

venn.diagram(summary_up_pathway, filename = 'up_pathway_fdr.png', imagetype = 'png',
             category.names = c("hypo_up" , "hypo_demeth" , "tcga_edgeR_up","tcga_wilcox_up","tcga_cor_neg"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), cat.cex = 0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), cex = 1, fontfamily = 'serif')


methy_delta_down_go_sig <- subset(methy_delta_down_go,subset = p.adjust<0.05)
hypo_down_go_sig <- subset(hypo_down_go,subset = p.adjust<0.05)
tcga_pos_go_sig <- subset(tcga_pos_go,subset = p.adjust<0.05)
tcga_edgeR_down_go_sig <-  subset(tcga_edgeR_down_go,subset = p.adjust<0.05)
tcga_wilcox_down_go_sig <-  subset(tcga_wilcox_down_go,subset = pvalue<0.05)

summary_down_pathway <- list(hypo_down_go_sig$Description,
                             methy_delta_down_go_sig$Description,
                             tcga_edgeR_down_go_sig$Description,
                             tcga_wilcox_down_go_sig$Description,
                             tcga_pos_go_sig$Description)
intersection_down_pathway <- Reduce(intersect,summary_down_pathway)
venn.diagram(summary_down_pathway, filename = 'down_pathway_fdr.png', imagetype = 'png',
             category.names = c("hypo_up" , "hypo_demeth" , "tcga_edgeR_up","tcga_wilcox_up","tcga_cor_neg"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), cat.cex = 0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#EAC7F2','#F75A6F'), cex = 1, fontfamily = 'serif')


metabolic_gene <- subset(tcga_edgeR,subset = V1 %in% genelist)

tcga_edgeR <- mutate(tcga_edgeR,matabolic_gene = case_when(V1 %in% genelist ~ 'metabolic_gene'))

tcga_edgeR$matabolic_gene <- replace(tcga_edgeR$matabolic_gene,is.na(tcga_edgeR$matabolic_gene),'other')


genelist <- c('CYP2A7/GLYAT/NR1I2/CYP1A1/CYP2A6/UGT1A4/GSTA2/GSTA3/CYP2A13/CYP1A2/CYP2G1P/CYP3A4/GSTA1/UGT1A9/CYP2B6/NAT2/CYP2D7/ACSM1/UGT2B28/ACSM2B/UGT1A1/CYP2F1/SULT1B1/CMBL/CYP2J2/AOC1/CES1/UGT2B11/GSTA5/CYP2C18/CES2/FMO3/CES3/UGT1A6')
genelist <- unlist(strsplit(genelist,"/"))
genelist2 <- c('OR56A3/OR8A1/OR5H8/OR5H2/OR8G5/OR8G3P/OR52E6/OR51B5/OR52N5/OR52E4/OR51Q1/OR51B4/GFY/OR8D1/OR2C3/OR10J5/OR6C70/OR1F1/OR6C2/OR52E5/OR52N1/OR51J1/OR5M11/OR4D10/OR5M8')
genelist2 <- unlist(strsplit(genelist2,"/"))

wilcox_OR <- subset(tcga_wilcox,subset=gene %in% genelist2)
wilcox_OR_exp_hypomean <- dplyr::select(wilcox_OR,gene,hypomean)
wilcox_OR_exp_hypermean <- dplyr::select(wilcox_OR,gene,hypermean)
colnames(wilcox_OR_exp_hypomean) <- c("gene","mean_cpm")
colnames(wilcox_OR_exp_hypermean) <- c("gene","mean_cpm")
wilcox_OR_exp_hypomean$group <- 'hypo'
wilcox_OR_exp_hypermean$group <- 'hyper'
wilcox_OR_sum <- rbind(wilcox_OR_exp_hypomean,wilcox_OR_exp_hypermean)
ggplot(wilcox_OR_sum,aes(x=group,y=mean_cpm,color=group))+
  geom_boxplot(width=0.5,outlier.size = 0)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("hypo","hyper")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =4 )+theme_bw() +geom_point(aes(fill=group, group=gene), size = 2,position = position_dodge(0.2), colour = "black")+
  geom_line(aes(group=gene), color="gray" ,position = position_dodge(0.2))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())








