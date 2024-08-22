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

venn.diagram(hypo_up_genes_list, filename = 'hypo_up.png', imagetype = 'png',
             category.names = c("hcc2" , "hcc8" , "hcc9"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex =0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')


hypo_down_genes_list <- list(hcc2_hypo_genes_down_sig$gene,hcc8_hypo_genes_down_sig$gene,hcc9_hypo_genes_down_sig$gene)
common_down_genelist <- Reduce(intersect,hypo_down_genes_list)
