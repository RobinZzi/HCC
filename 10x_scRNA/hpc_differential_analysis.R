
hpc_merge_10x <- subset(hpc_merge_harmony,subset = tech =="10x")




hcc3_10x_demeth$gene <- row.names(hcc3_10x_demeth) 
hcc28_10x_demeth$gene <- row.names(hcc28_10x_demeth) 
hcc29_10x_demeth$gene <- row.names(hcc29_10x_demeth) 



hcc3_10x_demeth_up <- subset(hcc3_10x_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc3_10x_demeth_down <- subset(hcc3_10x_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)

hcc28_10x_demeth_up <- subset(hcc28_10x_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc28_10x_demeth_down <- subset(hcc28_10x_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)

hcc29_10x_demeth_up <- subset(hcc29_10x_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc29_10x_demeth_down <- subset(hcc29_10x_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)




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

up_10x_venn <- list(hcc3_10x_demeth_up$gene,hcc28_10x_demeth_up$gene,hcc29_10x_demeth_up$gene)


venn.diagram(up_10x_venn, filename = 'up_10x.png', imagetype = 'png', ,category.names = c("hcc3" , "hcc28" , "hcc29"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

down_10x_venn <- list(hcc3_10x_demeth_down$gene,hcc28_10x_demeth_down$gene,hcc29_10x_demeth_down$gene)

venn.diagram(down_10x_venn, filename = 'down_10x.png', imagetype = 'png', ,category.names = c("hcc3" , "hcc28" , "hcc29"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')


hpc_merge_strt <- subset(hpc_merge_harmony,subset = tech =="strt")



hpc_merge_strt <- subset(hpc_merge_harmony,subset = tech =="strt")

hcc3_strt_demeth <- FindMarkers(hpc_merge_strt,ident.1 = "hcc3_demeth",ident.2 = "hcc3_admeth",group.by = "hpc_meth_type")
hcc28_strt_demeth <- FindMarkers(hpc_merge_strt,ident.1 = "hcc28_demeth",ident.2 = "hcc3_admeth",group.by = "hpc_meth_type")
hcc29_strt_demeth <- FindMarkers(hpc_merge_strt,ident.1 = "hcc29_demeth",ident.2 = "hcc3_admeth",group.by = "hpc_meth_type")



hcc3_strt_demeth$gene <- row.names(hcc3_strt_demeth) 
hcc28_strt_demeth$gene <- row.names(hcc28_strt_demeth) 
hcc29_strt_demeth$gene <- row.names(hcc29_strt_demeth) 



hcc3_strt_demeth_up <- subset(hcc3_strt_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc3_strt_demeth_down <- subset(hcc3_strt_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)

hcc28_strt_demeth_up <- subset(hcc28_strt_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc28_strt_demeth_down <- subset(hcc28_strt_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)

hcc29_strt_demeth_up <- subset(hcc29_strt_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc29_strt_demeth_down <- subset(hcc29_strt_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)




common_strt_down <- intersect(hcc3_strt_demeth_down$gene,hcc28_strt_demeth_down$gene)
common_strt_down <- intersect(common_strt_down ,hcc29_strt_demeth_down$gene)


common_strt_up <- intersect(hcc3_strt_demeth_up$gene,hcc28_strt_demeth_up$gene)
common_strt_up <- intersect(common_strt_up,hcc29_strt_demeth_up$gene)









hpc_common_strt_up_go <- enrichGO(gene  = common_strt_up,
                                  OrgDb      = org.Hs.eg.db,
                                  keyType    = 'SYMBOL',
                                  ont        = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
hpc_common_strt_up_go <- as.data.frame(hpc_common_strt_up_go@result)
hpc_common_strt_up_go [,"logp"] <- -log10(hpc_common_strt_up_go$pvalue)
ggplot(data = hpc_common_strt_up_go[1:20,])+
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








hpc_common_strt_down_go <- enrichGO(gene  = common_strt_down,
                                    OrgDb      = org.Hs.eg.db,
                                    keyType    = 'SYMBOL',
                                    ont        = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05)
hpc_common_strt_down_go <- as.data.frame(hpc_common_strt_down_go@result)
hpc_common_strt_down_go [,"logp"] <- -log10(hpc_common_strt_down_go$pvalue)
ggplot(data = hpc_common_strt_down_go[1:20,])+
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

up_strt_venn <- list(hcc3_strt_demeth_up$gene,hcc28_strt_demeth_up$gene,hcc29_strt_demeth_up$gene)


venn.diagram(up_strt_venn, filename = 'up_strt.png', imagetype = 'png', ,category.names = c("hcc3" , "hcc28" , "hcc29"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

down_strt_venn <- list(hcc3_strt_demeth_down$gene,hcc28_strt_demeth_down$gene,hcc29_strt_demeth_down$gene)

venn.diagram(down_strt_venn, filename = 'down_strt.png', imagetype = 'png', ,category.names = c("hcc3" , "hcc28" , "hcc29"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')







write.table(hcc3_10x_demeth_up,"hcc3_10x_demeth_up.txt")
write.table(hcc3_10x_demeth_down,"hcc3_10x_demeth_down.txt")
write.table(hcc28_10x_demeth_up,"hcc28_10x_demeth_up.txt")
write.table(hcc28_10x_demeth_down,"hcc28_10x_demeth_down.txt")
write.table(hcc29_10x_demeth_up,"hcc29_10x_demeth_up.txt")
write.table(hcc29_10x_demeth_down,"hcc29_10x_demeth_down.txt")









write.table(hcc3_strt_demeth_up,"hcc3_strt_demeth_up.txt")
write.table(hcc3_strt_demeth_down,"hcc3_strt_demeth_down.txt")
write.table(hcc28_strt_demeth_up,"hcc28_strt_demeth_up.txt")
write.table(hcc28_strt_demeth_down,"hcc28_strt_demeth_down.txt")
write.table(hcc29_strt_demeth_up,"hcc29_strt_demeth_up.txt")
write.table(hcc29_strt_demeth_down,"hcc29_strt_demeth_down.txt")


