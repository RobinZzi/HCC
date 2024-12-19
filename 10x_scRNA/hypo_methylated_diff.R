table(HCC2_HPC$methinfo)




DotPlot(demeth_hpc,features="SNHG6",group.by="patient_pt")
VlnPlot(demeth_hpc,features="SNHG6",group.by="patient_pt")
VlnPlot(HPC,features="SNHG6",group.by="patient_pt")+NoLegend()
VlnPlot(demeth_hpc,features="GADD45A",group.by="patient_pt")+NoLegend()


demeth_hpc_pt_info <- demeth_hpc$patient_pt

demeth_hpc_meth_info = case_when(
  demeth_hpc_pt_info  %in% c("HCC2_PT3","HCC2_NT")~"hcc2_normal_meth",
  demeth_hpc_pt_info  %in% c("HCC2_PT1","HCC2_PT2")~"hcc2_hypo_meth",
  demeth_hpc_pt_info  %in% c("HCC8_PT1","HCC8_PT2","HCC8_PT4")~"hcc8_hypo_meth",
  demeth_hpc_pt_info  %in% c("HCC9_PT1","HCC9_PT3","HCC9_PT4")~"hcc9_hypo_meth",
  demeth_hpc_pt_info  %in% c("HCC8_NT")~"hcc8_normal_meth",
  demeth_hpc_pt_info  %in% c("HCC9_NT")~"hcc9_normal_meth",
  TRUE ~ as.character(demeth_hpc_pt_info))

demeth_hpc$meth_info <- demeth_hpc_meth_info

DimPlot(demeth_hpc,group.by = "meth_info",label=T)


hcc2_hypo_genes <- FindMarkers(demeth_hpc,ident.1 = "hcc2_hypo_meth",ident.2 = "hcc2_normal_meth",group.by="meth_info")
hcc2_hypo_genes$gene <- row.names(hcc2_hypo_genes)
hcc2_hypo_genes_up_sig <- subset(hcc2_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC>0.5)
hcc2_hypo_genes_down_sig <- subset(hcc2_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC<0)

write.table(hcc2_hypo_genes_up_sig,"hcc2_hypo_up_genes.txt")


hcc8_hypo_genes <- FindMarkers(demeth_hpc,ident.1 = "hcc8_hypo_meth",ident.2 = "hcc8_normal_meth",group.by="meth_info")
hcc8_hypo_genes$gene <- row.names(hcc8_hypo_genes)
hcc8_hypo_genes_sig <- subset(hcc8_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC>0)
hcc8_hypo_genes_up_sig <- subset(hcc8_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC>0)
hcc8_hypo_genes_down_sig <- subset(hcc8_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC<0)

hcc9_hypo_genes <- FindMarkers(demeth_hpc,ident.1 = "hcc9_hypo_meth",ident.2 = "hcc9_normal_meth",group.by="meth_info")
hcc9_hypo_genes$gene <- row.names(hcc9_hypo_genes)
hcc9_hypo_genes_sig <- subset(hcc9_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC>0)
hcc9_hypo_genes_up_sig <- subset(hcc9_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC>0)
hcc9_hypo_genes_down_sig <- subset(hcc9_hypo_genes,subset = p_val_adj<0.05 & avg_log2FC<0)


hypo_up_genes_list <- list(hcc2_hypo_genes_up_sig$gene,hcc8_hypo_genes_up_sig$gene,hcc9_hypo_genes_up_sig$gene)
common_up_genelist <- Reduce(intersect,hypo_up_genes_list)


top_neg <- fread("cor_data_df_neg_top.txt")

hypo_up_genes_tcga_list <- list(hcc2_hypo_genes_up_sig$gene,hcc8_hypo_genes_up_sig$gene,
                                hcc9_hypo_genes_up_sig$gene,top_neg$symbol)


venn.diagram(hypo_up_genes_tcga_list, filename = 'hypo_up_0912_includetcga.png', imagetype = 'png',
             category.names = c("hcc2" , "hcc8" , "hcc9","tcga_neg"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#fc8b93'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#fc8b93'), cat.cex =1, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#fc8b93'), cex = 1.5, fontfamily = 'serif')



venn.diagram(hypo_up_genes_list, filename = 'hypo_up_1215.png', imagetype = 'png',
             category.names = c("hcc2" , "hcc8" , "hcc9"),
             fill = c('#ee9caa', '#13a1c1', '#472767'), alpha = 0.50, 
             cat.col = c('#ee9caa', '#13a1c1', '#472767'), cat.cex =0, cat.fontfamily = 'serif',
             col = c('#ee9caa', '#13a1c1', '#472767'), cex = 1.5, fontfamily = 'serif')


hypo_down_genes_list <- list(hcc2_hypo_genes_down_sig$gene,hcc8_hypo_genes_down_sig$gene,hcc9_hypo_genes_down_sig$gene)
common_down_genelist <- Reduce(intersect,hypo_down_genes_list)



tf_list <- list(hcc2_hypo_genes_up_sig$gene,hcc8_hypo_genes_up_sig$gene,
                                           hcc9_hypo_genes_up_sig$gene)
tf_result <- Reduce(intersect,tf_list)

tf_go <- enrichGO(gene  = tf_result,
                                OrgDb      = org.Hs.eg.db,
                                keyType    = 'SYMBOL',
                                ont        = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
tf_go <- as.data.frame(tf_go@result)
tf_go [,"LogFDR"] <- -log10(tf_go$p.adjust)
tf_go [,"logp"] <- -log10(tf_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
tf_enrichment_fold <- apply(tf_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
tf_go$EF <- tf_enrichment_fold

#Draw Bubble Plot
ggplot(data = tf_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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
       title="Up-Regulated Pathways in de-methylated samples of HCC")


ts_list <- list(hcc2_hypo_genes_up_sig$gene,top_neg$symbol)
ts_result <- Reduce(intersect,ts_list)

ts_go <- enrichGO(gene  = ts_result,
                  OrgDb      = org.Hs.eg.db,
                  keyType    = 'SYMBOL',
                  ont        = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
ts_go <- as.data.frame(ts_go@result)
ts_go [,"LogFDR"] <- -log10(ts_go$p.adjust)
ts_go [,"logp"] <- -log10(ts_go$pvalue)
#Use GeneRatio and BgRatio to calculate Enrichment_Fold(EF = GR/BR)
ts_enrichment_fold <- apply(ts_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
ts_go$EF <- ts_enrichment_fold

#Draw Bubble Plot
ggplot(data = ts_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
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
       title="Up-Regulated Pathways in de-methylated samples of HCC2 & TCGA")

