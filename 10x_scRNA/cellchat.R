
install.packages(c('NMF','circlize'))
BiocManager::install("ComplexHeatmap")
BiocManager::install("BiocNeighbors")
# 包的路径
devtools::install_local("/storage/zhangyanxiaoLab/zhangliwen/src/CellChat-master.zip")

library(CellChat)
CellChatDB <- CellChatDB.human

library(Seurat)

merge_cmn <- subset(bigseu, subset = cntype %in% "CMN")
merge_sn <- subset(bigseu, subset = cntype %in% "SN")
merge_pt <- subset(bigseu, subset = group %in% "P")
merge_st <- subset(bigseu, subset = group %in% "S")



cmn_data_input <- GetAssayData(merge_cmn, assay = "RNA", slot = "data")
cmn_identityentity <- subset(merge_cmn@meta.data, select = "new.ident")
cmn_cellchat <- createCellChat(object = cmn_data_input, meta = cmn_identity,  group.by = "new.ident")
cmn_CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cmn_cellchat@DB <- cmn_CellChatDB.use
cmn_cellchat <- subsetData(cmn_cellchat)
cmn_cellchat <- identifyOverExpressedGenes(cmn_cellchat)
cmn_cellchat <- identifyOverExpressedInteractions(cmn_cellchat)
cmn_cellchat <- projectData(cmn_cellchat, PPI.human)
cmn_cellchat <- computeCommunProb(cmn_cellchat, raw.use = TRUE)
cmn_cellchat <- filterCommunication(cmn_cellchat, min.cells = 3)
cmn_df.net <- subsetCommunication(cmn_cellchat)
cmn_cellchat <- computeCommunProbPathway(cmn_cellchat)
cmn_cellchat <- aggregateNet(cmn_cellchat)
cmn_groupSize <- as.numeric(table(cmn_cellchat@idents))


CairoPNG("CMN-cellchat.png",xpd=TRUE,height = 800,width = 800)
netVisual_circle(cmn_cellchat@net$count, vertex.weight = cmn_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")




sn_data_input <- GetAssayData(merge_sn, assay = "RNA", slot = "data")
sn_identity <- subset(merge_sn@meta.data, select = "new.ident")
sn_identity = droplevels(sn_identity, exclude = setdiff(levels(sn_identity),unique(sn_identity)))
sn_cellchat <- createCellChat(object = sn_data_input, meta = sn_identity,  group.by = "new.ident")
sn_CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
sn_cellchat@DB <- sn_CellChatDB.use
sn_cellchat <- subsetData(sn_cellchat)
sn_cellchat <- identifyOverExpressedGenes(sn_cellchat)
sn_cellchat <- identifyOverExpressedInteractions(sn_cellchat)
sn_cellchat <- projectData(sn_cellchat, PPI.human)
sn_cellchat <- computeCommunProb(sn_cellchat, raw.use = TRUE)
sn_cellchat <- filterCommunication(sn_cellchat, min.cells = 3)
sn_df.net <- subsetCommunication(sn_cellchat)
sn_cellchat <- computeCommunProbPathway(sn_cellchat)
sn_cellchat <- aggregateNet(sn_cellchat)
sn_groupSize <- as.numeric(table(sn_cellchat@idents))


CairoPNG("SN-cellchat.png",xpd=TRUE,height = 800,width = 800)
netVisual_circle(sn_cellchat@net$count, vertex.weight = sn_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")


pt_data_input <- GetAssayData(merge_pt, assay = "RNA", slot = "data")
pt_identity <- subset(merge_pt@meta.data, select = "new.ident")
pt_identity = droplevels(pt_identity, exclude = setdiff(levels(pt_identity),unique(pt_identity)))
pt_cellchat <- createCellChat(object = pt_data_input, meta = pt_identity,  group.by = "new.ident")
pt_CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
pt_cellchat@DB <- pt_CellChatDB.use
pt_cellchat <- subsetData(pt_cellchat)
pt_cellchat <- identifyOverExpressedGenes(pt_cellchat)
pt_cellchat <- identifyOverExpressedInteractions(pt_cellchat)
pt_cellchat <- projectData(pt_cellchat, PPI.human)
pt_cellchat <- computeCommunProb(pt_cellchat, raw.use = TRUE)
pt_cellchat <- filterCommunication(pt_cellchat, min.cells = 3)
pt_df.net <- subsetCommunication(pt_cellchat)
pt_cellchat <- computeCommunProbPathway(pt_cellchat)
pt_cellchat <- aggregateNet(pt_cellchat)
pt_groupSize <- as.numeric(table(pt_cellchat@idents))

CairoPNG("PT-cellchat.png",xpd=TRUE,height = 800,width = 800)
netVisual_circle(pt_cellchat@net$count, vertex.weight = pt_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")


st_data_input <- GetAssayData(merge_st, assay = "RNA", slot = "data")
st_identity <- subset(merge_st@meta.data, select = "new.ident")
st_cellchat <- createCellChat(object = st_data_input, meta = st_identity,  group.by = "new.ident")
st_CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
st_cellchat@DB <- st_CellChatDB.use
st_cellchat <- subsetData(st_cellchat)
st_cellchat <- identifyOverExpressedGenes(st_cellchat)
st_cellchat <- identifyOverExpressedInteractions(st_cellchat)
st_cellchat <- projectData(st_cellchat, PPI.human)
st_cellchat <- computeCommunProb(st_cellchat, raw.use = TRUE)
st_cellchat <- filterCommunication(st_cellchat, min.cells = 3)
st_df.net <- subsetCommunication(st_cellchat)
st_cellchat <- computeCommunProbPathway(st_cellchat)
st_cellchat <- aggregateNet(st_cellchat)
st_groupSize <- as.numeric(table(st_cellchat@idents))


CairoPNG("ST-cellchat.png",xpd=TRUE,height = 800,width = 800)
netVisual_circle(st_cellchat@net$count, vertex.weight = st_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")














merge_cmn_drop <- subset(bigseu, subset = lib.method == "dropseq")
drop_data_input <- GetAssayData(merge_cmn_drop, assay = "RNA", slot = "data")
drop_identity <- subset(merge_cmn_drop@meta.data, select = "new.ident")
drop_cellchat <- createCellChat(object = drop_data_input, meta = drop_identity,  group.by = "new.ident")
drop_CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
drop_cellchat@DB <- drop_CellChatDB.use
drop_cellchat <- subsetData(drop_cellchat)
drop_cellchat <- identifyOverExpressedGenes(drop_cellchat)
drop_cellchat <- identifyOverExpressedInteractions(drop_cellchat)
drop_cellchat <- projectData(drop_cellchat, PPI.human)
drop_cellchat <- computeCommunProb(drop_cellchat, raw.use = TRUE)
drop_cellchat <- filterCommunication(drop_cellchat, min.cells = 3)
drop_df.net <- subsetCommunication(drop_cellchat)
drop_cellchat <- computeCommunProbPathway(drop_cellchat)
drop_cellchat <- aggregateNet(drop_cellchat)
drop_groupSize <- as.numeric(table(drop_cellchat@idents))


CairoPNG("HCC6-cellchat.png",xpd=TRUE,height = 800,width = 800)
netVisual_circle(drop_cellchat@net$count, vertex.weight = drop_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

cmn_sn_drop_object_list <- list(cmn = drop_cellchat,sn = sn_cellchat)
cmn_sn_drop_cellchat <- mergeCellChat(cmn_sn_drop_object_list, add.names = names(cmn_sn_drop_object_list),cell.prefix = TRUE)

rankNet(cmn_sn_drop_cellchat, mode = "comparison", stacked = F,do.stat = FALSE)



cmn_sn_object_list <- list(cmn = cmn_cellchat,sn = sn_cellchat)
cmn_sn_cellchat <- mergeCellChat(cmn_sn_object_list, add.names = names(cmn_sn_object_list),cell.prefix = TRUE)
cmn_sn_gg1 <- compareInteractions(cmn_sn_cellchat, show.legend = F, group = c(1:3))
cmn_sn_gg2 <- compareInteractions(cmn_sn_cellchat, show.legend = F, group = c(1:3), measure = "weight")
cmn_sn_gg1 + cmn_sn_gg2
cmn_sn_hm1 <- netVisual_heatmap(cmn_sn_cellchat)
cmn_sn_hm2 <- netVisual_heatmap(cmn_sn_cellchat, measure = "weight")
cmn_sn_hm1 + cmn_sn_hm2

netVisual_heatmap(cmn_sn_cellchat)
netVisual_heatmap(cmn_sn_cellchat, measure = "weight")
netVisual_diffInteraction(cmn_sn_cellchat, weight.scale = T, measure = "count.merged", label.edge = T)

pt_st_object_list <- list(pt = pt_cellchat, st = st_cellchat)
pt_st_cellchat <- mergeCellChat(pt_st_object_list, add.names = names(pt_st_object_list))
pt_st_gg1 <- compareInteractions(pt_st_cellchat, show.legend = F, group = c(1:2))
pt_st_gg2 <- compareInteractions(pt_st_cellchat, show.legend = F, group = c(1:2), measure = "weight")
pt_st_gg1 + pt_st_gg2



cmn_samples_cellchat <- data.frame()
merge_cmn <- subset(bigseu, subset = cntype %in% "CMN")
merge_cmn2 <- subset(merge_cmn, subset = group != "NT")
cmn_samples <- unique(merge_cmn2@meta.data$patient_pt)
cmn_samples_cellchat_pwlist <- list()
# 对每个样品来源，计算每种细胞类型的比例
for(s in cmn_samples){
  cells_of_sample <- subset(merge_cmn,subset=patient_pt == s)
  counts_of_sample <- GetAssayData(cells_of_sample, assay = "RNA", slot = "data")
  identity <- subset(cells_of_sample@meta.data, select = "new.ident")
  identity = droplevels(identity, exclude = setdiff(levels(identity),unique(identity)))
  cellchat <- createCellChat(object = counts_of_sample, meta = identity,  group.by = "new.ident")
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- cellchat %>% 
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    projectData(PPI.human) %>%
    computeCommunProb(raw.use = TRUE) %>%
    filterCommunication(min.cells = 3) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  prop_df <- data.frame(
    'sample' = s, 
    'interations' = sum(cellchat@net[["count"]]),
    'weight' = sum(cellchat@net[["weight"]]),
    'type' = 'cmn'
  )
  p <- rankNet(cellchat, mode = "single", stacked = F,do.stat = FALSE)
  cmn_samples_cellchat_pwlist[[s]] <- p$data
  cmn_samples_cellchat <- rbind(cmn_samples_cellchat, prop_df)
}



sn_samples_cellchat <- data.frame()
merge_sn <- subset(bigseu, subset = cntype %in% "SN")
merge_sn2 <- subset(merge_sn, subset = group != "NT")
sn_samples <- unique(merge_sn2@meta.data$patient_pt)
sn_samples_cellchat_pwlist <- list()
for(s in sn_samples){
  cells_of_sample <- subset(merge_sn,subset=patient_pt == s)
  counts_of_sample <- GetAssayData(cells_of_sample, assay = "RNA", slot = "data")
  identity <- subset(cells_of_sample@meta.data, select = "new.ident")
  identity = droplevels(identity, exclude = setdiff(levels(identity),unique(identity)))
  cellchat <- createCellChat(object = counts_of_sample, meta = identity,  group.by = "new.ident")
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- cellchat %>% 
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    projectData(PPI.human) %>%
    computeCommunProb(raw.use = TRUE) %>%
    filterCommunication(min.cells = 3) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  prop_df <- data.frame(
    'sample' = s, 
    'interactions' = sum(cellchat@net[["count"]]),
    'weight' = sum(cellchat@net[["weight"]]),
    'type' = 'sn'
  )
  p <- rankNet(cellchat, mode = "single", stacked = F,do.stat = FALSE)
  sn_samples_cellchat_pwlist[[s]] <- p$data
  sn_samples_cellchat <- rbind(sn_samples_cellchat, prop_df)
}

pt_samples_cellchat <- data.frame()
merge_pt <- subset(bigseu, subset = group %in% "P")
pt_samples <- unique(merge_pt@meta.data$patient_pt)

for(s in pt_samples){
  cells_of_sample <- subset(merge_pt,subset=patient_pt == s)
  counts_of_sample <- GetAssayData(cells_of_sample, assay = "RNA", slot = "data")
  identity <- subset(cells_of_sample@meta.data, select = "new.ident")
  identity = droplevels(identity, exclude = setdiff(levels(identity),unique(identity)))
  cellchat <- createCellChat(object = counts_of_sample, meta = identity,  group.by = "new.ident")
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- cellchat %>% 
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    projectData(PPI.human) %>%
    computeCommunProb(raw.use = TRUE) %>%
    filterCommunication(min.cells = 3) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  prop_df <- data.frame(
    'sample' = s, 
    'interactions' = sum(cellchat@net[["count"]]),
    'weight' = sum(cellchat@net[["weight"]]),
    'type' = 'pt'
  )
  pt_samples_cellchat <- rbind(pt_samples_cellchat, prop_df)
}

st_samples_cellchat <- data.frame()
merge_st <- subset(bigseu, subset = group %in% "S")
st_samples <- unique(merge_st@meta.data$patient_pt)

for(s in st_samples){
  cells_of_sample <- subset(merge_st,subset=patient_pt==s)
  counts_of_sample <- GetAssayData(cells_of_sample, assay = "RNA", slot = "data")
  identity <- subset(cells_of_sample@meta.data, select = "new.ident")
  identity = droplevels(identity, exclude = setdiff(levels(identity),unique(identity)))
  cellchat <- createCellChat(object = counts_of_sample, meta = identity,  group.by = "new.ident")
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- cellchat %>% 
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    projectData(PPI.human) %>%
    computeCommunProb(raw.use = TRUE) %>%
    filterCommunication(min.cells = 3) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  prop_df <- data.frame(
    'sample' = s, 
    'interactions' = sum(cellchat@net[["count"]]),
    'weight' = sum(cellchat@net[["weight"]]),
    'type' = 'st'
  )
  st_samples_cellchat <- rbind(st_samples_cellchat, prop_df)
}

dropseq_list <- unique(subset(bigseu,subset=lib.method=='dropseq')$patient_pt)




cmnsn_cellchat_sum <- rbind(cmn_samples_cellchat,sn_samples_cellchat)
ptst_cellchat_sum <- rbind(pt_samples_cellchat,st_samples_cellchat)


cmnsn_cellchat_sum_drop <- subset(cmnsn_cellchat_sum,subset=sample%in%dropseq_list)

ggplot(cmnsn_cellchat_sum,aes(x=type,y=interactions,color=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 20)+
  labs(x = "type", 
       y= "interaction",
       title="10x excluded")+
  scale_color_manual(values=c("sn"="#3f72af","cmn"="#d72323"))+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y =1200 )+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

ggplot(cmnsn_cellchat_sum_drop,aes(x=type,y=interactions,color=type))+
  geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.1,shape = 20)+
  labs(x = "type", 
       y= "interaction",
       title="drop-seq only")+
  scale_color_manual(values=c("sn"="#3f72af","cmn"="#d72323"))+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",
                     label.y =600 )+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())


ggplot(cmnsn_cellchat_sum,aes(x=type,y=interactions,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="cell_chat")+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y =80 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",color = 'black',size=12),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


ggplot(cmnsn_cellchat_sum_drop,aes(x=type,y=interactions,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="cell_chat")+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y =700 )+theme_bw()

ggplot(cmnsn_cellchat_sum_drop,aes(x=type,y=interactions,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="cell_chat")+
  scale_fill_manual(values=c("sn"="#3f72af","cmn"="#d72323"))+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y =600 )+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

ggplot(cmnsn_cellchat_sum_drop,aes(x=type,y=weight,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="cell_chat")+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y =20 )+theme_bw()



ggplot(ptst_cellchat_sum,aes(x=type,y=interactions,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="cell_chat")+
  stat_compare_means(comparisons = list(c("pt","st")),
                     method = "t.test",label = "p.signif",
                     label.y =1300 )+theme_bw()


ggplot(ptst_cellchat_sum,aes(x=type,y=weight,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="cell_chat")+
  stat_compare_means(comparisons = list(c("pt","st")),
                     method = "t.test",label = "p.signif",
                     label.y =80)+theme_bw()

ggplot(cmn_samples_cellchat,aes(x=type,y=interations))+
  geom_boxplot()+geom_jitter()+
  labs(x = "type", 
       y= "interaction",
       title="cell_chat")

ggplot(cmn_samples_cellchat,aes(x=type,y=weight))+
  geom_boxplot()+geom_jitter()+
  labs(x = "type", 
       y= "weight",
       title="cell_chat")









CairoPNG("HCC3_PT2-cellchat.png",xpd=TRUE,height = 800,width = 800)
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count,
                 color.use = cellchat_color_sn,
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

cellchat_color_sn <- c("#E95C59","#57C3F3","#E59CC4","#53A85F", "#D6E7A3","#F1BB72","#E0D4CA",
                       "#23452F","#476D87","#AB3282","#BD956A")

cellchat_color_cmn <- c("#E5D2DD","#E95C59","#57C3F3","#585658","#E59CC4","#53A85F", "#D6E7A3","#F1BB72","#E0D4CA",
                       "#23452F","#476D87","#AB3282","#BD956A","#F3B1A0")


netVisual_circle(cellchat@net$count,
                 color.use = cellchat_color_cmn,
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")



c( "HPC"="#E0D4CA", 
   "Fibroblast"="#9FA3A8", 
   "Endothelial cell" = "#585658",
   "Neutrophil" = "#8C549C",
   "PlasmaB cell" = "#BD956A", 
   "B cell" = "#23452F",
   "Proliferative T" = "#AB3282",
   "CD8+ memory" = "#E59CC4",
   "CD8+ exhausted" = "#E95C59",
   "CD8+ cytotoxic" = "#476D87",
   "CD4+ memory" = "#57C3F3",
   "CD4+ Treg" = "#D6E7A3",
   "Mast cell" = "#F3B1A0",
   "NK"="#F1BB72",
   "Dendritic cell"="#53A85F",
   "Macrophage"="#E5D2DD")

