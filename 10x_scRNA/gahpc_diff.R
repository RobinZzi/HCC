library(harmony)

hpc_merge <- readRDS("~/projects/hcc/analysis/merged_scrna/hpc_merge.Rds")


hpc_merge <- SCTransform(hpc_merge,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
hpc_merge <- RunPCA(hpc_merge, verbose = FALSE)
hpc_merge <- RunUMAP(hpc_merge, dims = 1:30, verbose = FALSE)
DimPlot(hpc_merge,group.by = "patient")
DimPlot(hpc_merge,group.by = "lib.method")
hpc_merge <- RunHarmony(hpc_merge, "lib.method")
hpc_merge <- RunPCA(hpc_merge, verbose = FALSE)
hpc_merge <- RunUMAP(hpc_merge,  dims = 1:15,reduction = "harmony")
hpc_merge <- FindNeighbors(hpc_merge, dims = 1:30, verbose = FALSE)
hpc_merge <- FindClusters(hpc_merge, verbose = FALSE)
DimPlot(hpc_merge,group.by = "lib.method")
DimPlot(hpc_merge,group.by = "patient")
DimPlot(hpc_merge,group.by = "subtype_info")
DimPlot(hpc_merge,group.by = "group")


hpc_merge_t <-subset(hpc_merge,subset=group%in% c("P","S"))


ntpt_diff <- FindMarkers(hpc_merge,ident.1 = "NT",group.by = "group" )


ntpt_diff <- FindMarkers(hpc_merge,ident.1 = "NT",group.by = "group" )

ntpt_diff_sig <- subset(ntpt_diff,subset = p_val_adj < 0.05)
ntpt_diff_sig_up <- subset(ntpt_diff_sig,subset = avg_log2FC > 1)
ntpt_diff_sig_down <- subset(ntpt_diff_sig,subset = avg_log2FC < -1)
ntpt_diff <- mutate(ntpt_diff,state = case_when(avg_log2FC>1&p_val_adj < 0.05 ~ "Up-regulated", # 上调
                                                avg_log2FC< -1&p_val_adj < 0.05 ~ "Down-regulated", # 下调
                                                TRUE ~ "Unchanged"))


ntpt_up_go <- enrichGO(gene  = row.names(ntpt_diff_sig_up),
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
ntpt_up_go <- as.data.frame(ntpt_up_go@result)
ntpt_up_go [,"logp"] <- -log10(ntpt_up_go$pvalue)
ntpt_up_go[,11:12] <- as.numeric(str_split_fixed(ntpt_up_go$GeneRatio,"/",2))
ntpt_up_go$GeneRatio <- ntpt_up_go[,11]/ntpt_up_go[,12]
ggplot(data = ntpt_up_go[1:10,],aes(y=reorder(Description,logp),x=logp))+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="steelblue",high="steelblue")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="pt_down")


ntpt_down_go <- enrichGO(gene  = row.names(ntpt_diff_sig_down),
                         OrgDb      = org.Hs.eg.db,
                         keyType    = 'SYMBOL',
                         ont        = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
ntpt_down_go <- as.data.frame(ntpt_down_go@result)
ntpt_down_go [,"logp"] <- -log10(ntpt_down_go$pvalue)
ntpt_down_go[,11:12] <- as.numeric(str_split_fixed(ntpt_down_go$GeneRatio,"/",2))
ntpt_down_go$GeneRatio <- ntpt_down_go[,11]/ntpt_down_go[,12]
ggplot(data = ntpt_down_go[1:10,],aes(y=reorder(Description,logp),x=logp))+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="red",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="pt_up")

ntpt_diff2 <- ntpt_diff
ntpt_diff2$avg_log2FC <- -ntpt_diff2$avg_log2FC
ggplot(ntpt_diff2, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(size = 0.4, aes(color = state)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("red","grey", "steelblue"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="pt vs nt diff gene")





cmnsn_diff <- FindMarkers(hpc_merge_t,ident.1 = "CMN",ident.2 = "SN",group.by = "subtype_info" )



cmnsn_diff <- FindMarkers(hpc_merge_t,ident.1 = "CMN",ident.2 = "SN",group.by = "subtype_info" )

cmnsn_diff_sig <- subset(cmnsn_diff,subset = p_val_adj < 0.05)
cmnsn_diff_sig_up <- subset(cmnsn_diff_sig,subset = avg_log2FC > 1)
cmnsn_diff_sig_down <- subset(cmnsn_diff_sig,subset = avg_log2FC < -1)
cmnsn_diff <- mutate(cmnsn_diff,state = case_when(avg_log2FC>1&p_val_adj < 0.05 ~ "Up-regulated", # 上调
                                                  avg_log2FC< -1&p_val_adj < 0.05 ~ "Down-regulated", # 下调
                                                  TRUE ~ "Unchanged"))


cmnsn_up_go <- enrichGO(gene  = row.names(cmnsn_diff_sig_up),
                        OrgDb      = org.Hs.eg.db,
                        keyType    = 'SYMBOL',
                        ont        = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
cmnsn_up_go <- as.data.frame(cmnsn_up_go@result)
cmnsn_up_go [,"logp"] <- -log10(cmnsn_up_go$pvalue)
cmnsn_up_go[,11:12] <- as.numeric(str_split_fixed(cmnsn_up_go$GeneRatio,"/",2))
cmnsn_up_go$GeneRatio <- cmnsn_up_go[,11]/cmnsn_up_go[,12]
ggplot(data = cmnsn_up_go[1:10,],aes(y=reorder(Description,logp),x=logp))+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="steelblue",high="steelblue")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="cmn_up")

cmnsn_down_go <- enrichGO(gene  = row.names(cmnsn_diff_sig_down),
                          OrgDb      = org.Hs.eg.db,
                          keyType    = 'SYMBOL',
                          ont        = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
cmnsn_down_go <- as.data.frame(cmnsn_down_go@result)
cmnsn_down_go [,"logp"] <- -log10(cmnsn_down_go$pvalue)
cmnsn_down_go[,11:12] <- as.numeric(str_split_fixed(cmnsn_down_go$GeneRatio,"/",2))
cmnsn_down_go$GeneRatio <- cmnsn_down_go[,11]/cmnsn_down_go[,12]
ggplot(data = cmnsn_down_go[1:10,],aes(y=reorder(Description,logp),x=logp))+
  geom_bar(aes(y=reorder(Description,logp),x=logp,fill=Count),stat='identity')+
  scale_fill_gradient(expression(Count),low="steelblue",high="steelblue")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="sn_up")


ggplot(cmnsn_diff, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(size = 0.4, aes(color = state)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue","grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="cmn vs sn diff gene")


save.image("ga_diff.Rdata")



