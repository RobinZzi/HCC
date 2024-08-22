hcc1_pt_markers <- FindMarkers(HCC1_HPC,ident.1 = "primary",ident.2 = "nt",group.by = "satellite_region" )
hcc1_pt_markers_sig <- subset(hcc1_pt_markers,subset = p_val_adj < 0.05)
hcc1_pt_markers_sig_up <- subset(hcc1_pt_markers_sig,subset = avg_log2FC > 0)
hcc1_pt_markers_sig_down <- subset(hcc1_pt_markers_sig,subset = avg_log2FC < 0)

hcc1_pt_up_go <- enrichGO(gene  = row.names(hcc1_pt_markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
hcc1_pt_up_go <- as.data.frame(hcc1_pt_up_go@result)
hcc1_pt_up_go [,"logp"] <- -log10(hcc1_pt_up_go$pvalue)
ggplot(data = hcc1_pt_up_go[1:20,])+
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
       title="hcc1_pt_up")

hcc1_pt_down_go <- enrichGO(gene  = row.names(hcc1_pt_markers_sig_down),
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
hcc1_pt_down_go <- as.data.frame(hcc1_pt_down_go@result)
hcc1_pt_down_go [,"logp"] <- -log10(hcc1_pt_down_go$pvalue)
ggplot(data = hcc1_pt_down_go[1:20,])+
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
       title="hcc1_pt_down")




hcc2_pt_markers <- FindMarkers(HCC2_HPC,ident.1 = "NT",group.by = "orig.ident" )
hcc2_pt_markers_sig <- subset(hcc2_pt_markers,subset = p_val_adj < 0.05)
hcc2_pt_markers_sig_up <- subset(hcc2_pt_markers_sig,subset = avg_log2FC < 0)
hcc2_pt_markers_sig_down <- subset(hcc2_pt_markers_sig,subset = avg_log2FC > 0)

hcc2_pt_up_go <- enrichGO(gene  = row.names(hcc2_pt_markers_sig_up),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
hcc2_pt_up_go <- as.data.frame(hcc2_pt_up_go@result)
hcc2_pt_up_go [,"logp"] <- -log10(hcc2_pt_up_go$pvalue)
ggplot(data = hcc2_pt_up_go[1:20,])+
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
       title="hcc2_pt_up")

hcc2_pt_down_go <- enrichGO(gene  = row.names(hcc2_pt_markers_sig_down),
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
hcc2_pt_down_go <- as.data.frame(hcc2_pt_down_go@result)
hcc2_pt_down_go [,"logp"] <- -log10(hcc2_pt_down_go$pvalue)
ggplot(data = hcc2_pt_down_go[1:20,])+
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
       title="hcc2_pt_down")


hcc3_pt_markers <- FindMarkers(HCC3_HPC,ident.1 = "primary",ident.2 = "nt",group.by = "satellite_region" )
hcc3_pt_markers_sig <- subset(hcc3_pt_markers,subset = p_val_adj < 0.05)
hcc3_pt_markers_sig_up <- subset(hcc3_pt_markers_sig,subset = avg_log2FC > 0)
hcc3_pt_markers_sig_down <- subset(hcc3_pt_markers_sig,subset = avg_log2FC < 0)

hcc3_pt_up_go <- enrichGO(gene  = row.names(hcc3_pt_markers_sig_up),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
hcc3_pt_up_go <- as.data.frame(hcc3_pt_up_go@result)
hcc3_pt_up_go [,"logp"] <- -log10(hcc3_pt_up_go$pvalue)
ggplot(data = hcc3_pt_up_go[1:20,])+
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
       title="hcc3_pt_up")

hcc3_pt_down_go <- enrichGO(gene  = row.names(hcc3_pt_markers_sig_down),
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
hcc3_pt_down_go <- as.data.frame(hcc3_pt_down_go@result)
hcc3_pt_down_go [,"logp"] <- -log10(hcc3_pt_down_go$pvalue)
ggplot(data = hcc3_pt_down_go[1:20,])+
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
       title="hcc3_pt_down")

hcc5_pt_markers <- FindMarkers(HCC5_HPC,ident.1 = "primary",ident.2 = "nt",group.by = "satellite_region" )
hcc5_pt_markers_sig <- subset(hcc5_pt_markers,subset = p_val_adj < 0.05)
hcc5_pt_markers_sig_up <- subset(hcc5_pt_markers_sig,subset = avg_log2FC > 0)
hcc5_pt_markers_sig_down <- subset(hcc5_pt_markers_sig,subset = avg_log2FC < 0)

hcc5_pt_up_go <- enrichGO(gene  = row.names(hcc5_pt_markers_sig_up),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
hcc5_pt_up_go <- as.data.frame(hcc5_pt_up_go@result)
hcc5_pt_up_go [,"logp"] <- -log10(hcc5_pt_up_go$pvalue)
ggplot(data = hcc5_pt_up_go[1:20,])+
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
       title="hcc5_pt_up")

hcc5_pt_down_go <- enrichGO(gene  = row.names(hcc5_pt_markers_sig_down),
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
hcc5_pt_down_go <- as.data.frame(hcc5_pt_down_go@result)
hcc5_pt_down_go [,"logp"] <- -log10(hcc5_pt_down_go$pvalue)
ggplot(data = hcc5_pt_down_go[1:20,])+
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
       title="hcc5_pt_down")



hcc6_pt_markers <- FindMarkers(HCC6_HPC,ident.1 = "primary",ident.2 = "nt",group.by = "satellite_region" )
hcc6_pt_markers_sig <- subset(hcc6_pt_markers,subset = p_val_adj < 0.05)
hcc6_pt_markers_sig_up <- subset(hcc6_pt_markers_sig,subset = avg_log2FC > 0)
hcc6_pt_markers_sig_down <- subset(hcc6_pt_markers_sig,subset = avg_log2FC < 0)

hcc6_pt_up_go <- enrichGO(gene  = row.names(hcc6_pt_markers_sig_up),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
hcc6_pt_up_go <- as.data.frame(hcc6_pt_up_go@result)
hcc6_pt_up_go [,"logp"] <- -log10(hcc6_pt_up_go$pvalue)
ggplot(data = hcc6_pt_up_go[1:20,])+
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
       title="hcc6_pt_up")

hcc6_pt_down_go <- enrichGO(gene  = row.names(hcc6_pt_markers_sig_down),
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
hcc6_pt_down_go <- as.data.frame(hcc6_pt_down_go@result)
hcc6_pt_down_go [,"logp"] <- -log10(hcc6_pt_down_go$pvalue)
ggplot(data = hcc6_pt_down_go[1:20,])+
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
       title="hcc6_pt_down")



hcc7_pt_markers <- FindMarkers(HCC7_HPC,ident.1 = "primary",ident.2 = "nt",group.by = "satellite_region" )
hcc7_pt_markers_sig <- subset(hcc7_pt_markers,subset = p_val_adj < 0.05)
hcc7_pt_markers_sig_up <- subset(hcc7_pt_markers_sig,subset = avg_log2FC > 0)
hcc7_pt_markers_sig_down <- subset(hcc7_pt_markers_sig,subset = avg_log2FC < 0)

hcc7_pt_up_go <- enrichGO(gene  = row.names(hcc7_pt_markers_sig_up),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
hcc7_pt_up_go <- as.data.frame(hcc7_pt_up_go@result)
hcc7_pt_up_go [,"logp"] <- -log10(hcc7_pt_up_go$pvalue)
ggplot(data = hcc7_pt_up_go[1:20,])+
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
       title="hcc7_pt_up")

hcc7_pt_down_go <- enrichGO(gene  = row.names(hcc7_pt_markers_sig_down),
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
hcc7_pt_down_go <- as.data.frame(hcc7_pt_down_go@result)
hcc7_pt_down_go [,"logp"] <- -log10(hcc7_pt_down_go$pvalue)
ggplot(data = hcc7_pt_down_go[1:20,])+
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
       title="hcc7_pt_down")


hcc8_pt_markers <- FindMarkers(HCC8_HPC,ident.1 = "NT",group.by = "orig.ident" )
hcc8_pt_markers_sig <- subset(hcc8_pt_markers,subset = p_val_adj < 0.05)
hcc8_pt_markers_sig_up <- subset(hcc8_pt_markers_sig,subset = avg_log2FC < 0)
hcc8_pt_markers_sig_down <- subset(hcc8_pt_markers_sig,subset = avg_log2FC > 0)
hcc8_pt_markers_sig$gene <- row.names(hcc8_pt_markers_sig)
hcc9_pt_markers_sig$gene <- row.names(hcc9_pt_markers_sig)


hcc8_pt_up_go <- enrichGO(gene  = row.names(hcc8_pt_markers_sig_up),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
hcc8_pt_up_go <- as.data.frame(hcc8_pt_up_go@result)
hcc8_pt_up_go [,"logp"] <- -log10(hcc8_pt_up_go$pvalue)
ggplot(data = hcc8_pt_up_go[1:20,])+
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
       title="hcc8_pt_up")

hcc8_pt_down_go <- enrichGO(gene  = row.names(hcc8_pt_markers_sig_down),
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
hcc8_pt_down_go <- as.data.frame(hcc8_pt_down_go@result)
hcc8_pt_down_go [,"logp"] <- -log10(hcc8_pt_down_go$pvalue)
ggplot(data = hcc8_pt_down_go[1:20,])+
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
       title="hcc8_pt_down")

hcc9_pt_markers <- FindMarkers(HCC9_HPC,ident.1 = "NT",group.by = "orig.ident" )
hcc9_pt_markers_sig <- subset(hcc9_pt_markers,subset = p_val_adj < 0.05)
hcc9_pt_markers_sig_up <- subset(hcc9_pt_markers_sig,subset = avg_log2FC < 0)
hcc9_pt_markers_sig_down <- subset(hcc9_pt_markers_sig,subset = avg_log2FC > 0)

hcc9_pt_up_go <- enrichGO(gene  = row.names(hcc9_pt_markers_sig_up),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
hcc9_pt_up_go <- as.data.frame(hcc9_pt_up_go@result)
hcc9_pt_up_go [,"logp"] <- -log10(hcc9_pt_up_go$pvalue)
ggplot(data = hcc9_pt_up_go[1:20,])+
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
       title="hcc9_pt_up")

hcc9_pt_down_go <- enrichGO(gene  = row.names(hcc9_pt_markers_sig_down),
                            OrgDb      = org.Hs.eg.db,
                            keyType    = 'SYMBOL',
                            ont        = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
hcc9_pt_down_go <- as.data.frame(hcc9_pt_down_go@result)
hcc9_pt_down_go [,"logp"] <- -log10(hcc9_pt_down_go$pvalue)
ggplot(data = hcc9_pt_down_go[1:20,])+
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
       title="hcc9_pt_down")


ptvsnt_up_sum <- Reduce(union,list(row.names(hcc1_pt_markers_sig_up),row.names(hcc2_pt_markers_sig_up),row.names(hcc3_pt_markers_sig_up),
                           row.names(hcc5_pt_markers_sig_up),row.names(hcc6_pt_markers_sig_up),row.names(hcc7_pt_markers_sig_up),
                           row.names(hcc8_pt_markers_sig_up),row.names(hcc9_pt_markers_sig_up)))

ptvsnt_up_sum_list <- list(hcc1 = row.names(hcc1_pt_markers_sig_up),hcc2 =row.names(hcc2_pt_markers_sig_up),hcc3 =row.names(hcc3_pt_markers_sig_up),
                           hcc5 =row.names(hcc5_pt_markers_sig_up),hcc6 =row.names(hcc6_pt_markers_sig_up),hcc7 =row.names(hcc7_pt_markers_sig_up),
                           hcc8 =row.names(hcc8_pt_markers_sig_up),hcc9 =row.names(hcc9_pt_markers_sig_up))

VlnPlot(HPC,features = c("PKHD1","SYNE2","NME5","DCDC2","ONECUT1","DYNC2H1","WWTR1","CFAP221","LCA5","RAB3IP","DNAAF4"),group.by = "group",split.by = "patient")

ptvsnt_down_sum <- Reduce(union,list(row.names(hcc1_pt_markers_sig_down),row.names(hcc2_pt_markers_sig_down),row.names(hcc3_pt_markers_sig_down),
                           row.names(hcc5_pt_markers_sig_down),row.names(hcc6_pt_markers_sig_down),row.names(hcc7_pt_markers_sig_down),
                           row.names(hcc8_pt_markers_sig_down),row.names(hcc9_pt_markers_sig_down)))


hpc_mtx <- GetAssayData(object = HPC, slot = "data")

hpc_diff_genes <- union(row.names(hpc_markers_sig_up),row.names(hpc_markers_sig_down))

hpc_mtx_up <- hpc_mtx[row.names(hcc1_pt_markers_sig_up[1:20,]),]
hpc_mtx_down <- hpc_mtx[row.names(hcc1_pt_markers_sig_down[1:20,]),]

hpc_genelist <- union(row.names(hpc_markers_sig_up[1:20,]),row.names(hpc_markers_sig_down[1:20,]))
hpc_mtx_diff <- hpc_mtx[hpc_genelist,]
hpc_mtx_colanno <- as.data.frame(HPC@meta.data$group)
row.names(hpc_mtx_colanno) <- row.names(HPC@meta.data)
colnames(hpc_mtx_colanno) <-"group"
hpc_mtx_colanno <- arrange(hpc_mtx_colanno, group)
pheatmap(hpc_mtx_diff[,row.names(hpc_mtx_colanno)],show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,
         annotation_col = hpc_mtx_colanno)


pheatmap(hpc_mtx_up)
DoHeatmap(HPC, slot = "data",features = row.names(hcc1_pt_markers_sig_up[1:20,]),group.by ="group" )
DoHeatmap(HPC, slot = "data",features = row.names(hcc9_pt_markers_sig_down[1:20,]),group.by ="group" )
