rm(list = ls())
library(data.table)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("~/projects/hcc/analysis/meth_diff_gene")
save.image("meth_diff_gene.Rdata")


hcc1_uplist <- fread("pt1_pt2_up.txt")
hcc1_downlist <- fread("pt1_pt2_down.txt")


TCGA_pos_re <- fread("cor_data_df_pos_top.txt")
TCGA_neg_re <- fread("cor_data_df_neg_top.txt")

TCGA_diff <- fread("tcga_meth-type-based_wilcox-test_diff-gene_overlapwithhcc_result.txt")
TCGA_diff_hypo_up <- subset(TCGA_diff,subset = log2FC > 0)
TCGA_diff_hypo_down <- subset(TCGA_diff,subset = log2FC < 0)



intersect(intersect(TCGA_pos_re$symbol,hcc1_downlist$x),TCGA_diff_hypo_down$gene)
intersect(intersect(TCGA_neg_re$symbol,hcc1_uplist$x),TCGA_diff_hypo_up$gene)

up_venn <- list(TCGA_neg_re$symbol,hcc1_uplist$x,TCGA_diff_hypo_up$gene)
down_venn <- list(TCGA_pos_re$symbol,hcc1_downlist$x,TCGA_diff_hypo_down$gene)

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
       title="Up-Regulated Pathways in de-methylated samples")

ggplot(data = hypo_up_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
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
       title="Up-Regulated Pathways in sh-SNHG6 Huh7")





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
       title="Down-Regulated Pathways in de-methylated samples")



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


hcc3_10x_demeth_down <- fread("hcc3_demeth_down.txt")
hcc28_10x_demeth_down <- fread("hcc28_demeth_down.txt")
hcc29_10x_demeth_down <- fread("hcc29_demeth_down.txt")
hcc3_10x_demeth_up <- fread("hcc3_demeth_up.txt")
hcc28_10x_demeth_up <- fread("hcc28_demeth_up.txt")
hcc29_10x_demeth_up <- fread("hcc29_demeth_up.txt")


common_10x_down <- intersect(hcc3_10x_demeth_down$gene,hcc28_10x_demeth_down$gene)
common_10x_down <- intersect(common_10x_down ,hcc29_10x_demeth_down$gene)

common_10x_up <- intersect(hcc3_10x_demeth_up$gene,hcc28_10x_demeth_up$gene)
common_10x_up <- intersect(common_10x_up,hcc29_10x_demeth_up$gene)



hpc_common_10x_up_go <- enrichGO(gene  = common_10x_up,
                                 OrgDb      = org.Hs.eg.db,
                                 keyType    = 'SYMBOL',
                                 ont        = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
hpc_common_10x_up_go <- as.data.frame(hpc_common_10x_up_go@result)
hpc_common_10x_up_go [,"logp"] <- -log10(hpc_common_10x_up_go$pvalue)
ggplot(data = hpc_common_10x_up_go[1:20,])+
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
       title="Up Pathways in demeth")





hpc_common_10x_down_go <- enrichGO(gene  = common_10x_down,
                                   OrgDb      = org.Hs.eg.db,
                                   keyType    = 'SYMBOL',
                                   ont        = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05)
hpc_common_10x_down_go <- as.data.frame(hpc_common_10x_down_go@result)
hpc_common_10x_down_go [,"logp"] <- -log10(hpc_common_10x_down_go$pvalue)
ggplot(data = hpc_common_10x_down_go[1:20,])+
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
       title="Down Pathways in demeth")             



methy_delta <- fread("methy_delta_screen.txt")


methy_delta_up <- subset(methy_delta,subset=RNAlog2FC>0&DNA_meth_dt<0&FDR<0.05)
methy_delta_down <- subset(methy_delta,subset=RNAlog2FC<0&DNA_meth_dt>0&FDR<0.05)


trio_up <- intersect(methy_delta_up$gene,hcc1_uplist$x)
tcga_up <- intersect(TCGA_neg_re$symbol,TCGA_diff_hypo_up$gene)



common_up <- intersect(trio_up,tcga_up)









trio_up_go <- enrichGO(gene  = trio_up,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
trio_up_go <- as.data.frame(trio_up_go@result)
trio_up_go [,"LogFDR"] <- -log10(trio_up_go$p.adjust)
trio_up_go [,"logp"] <- -log10(trio_up_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
trio_up_enrichment_fold <- apply(trio_up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
trio_up_go$EF <- trio_up_enrichment_fold

#Draw Bubble Plot
ggplot(data = trio_up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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
       title="Up-Regulated Pathways in de-methylated samples")
trio_up_go_sig <- subset(trio_up_go,pvalue < 0.05)












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
       title="Up-Regulated Pathways in de-methylated samples")
tcga_up_go_sig <- subset(tcga_up_go,pvalue < 0.05)



common_pathway <- intersect(tcga_up_go_sig$Description,trio_up_go_sig$Description)






summary_list <- list(methy_delta_up$gene,hcc1_uplist$x,TCGA_neg_re$symbol,TCGA_diff_hypo_up$gene)
venn.diagram(summary_list, filename = 'sum_up.png', imagetype = 'png',category.names = c("hcc2_demethy" , "hcc2_rna_up" , "tcga_rna_up","tcga_methy_neg"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cat.cex = 0.8, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cex = 1.5, fontfamily = 'serif')




summary_down_list <- list(methy_delta_down$gene,hcc1_downlist$x,TCGA_pos_re$symbol,TCGA_diff_hypo_down$gene)
venn.diagram(summary_down_list, filename = 'sum_down.png', imagetype = 'png', ,category.names = c("hcc2_admethy" , "hcc2_rna_down" , "tcga_rna_down","tcga_methy_pos"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cat.cex = 0.8, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#87D1AC'), cex = 1.5, fontfamily = 'serif')
