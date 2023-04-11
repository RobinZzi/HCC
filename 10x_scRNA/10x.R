setwd("~/projects/hcc/analysis/10x_scRNA")
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(DoubletFinder)

rm(list=ls())





save.image("10x.Rdata")

hcc29.seu2 <- CreateSeuratObject(counts = hcc29_cnt,
                                meta.data = hcc29_meta)
hcc29.seu2 <- PercentageFeatureSet(hcc29.seu2, pattern = "^MT-", col.name = "percent.mt")
hcc29.seu2 <- subset(hcc29.seu2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hcc29.seu <- NormalizeData(hcc29.seu, normalization.method = "LogNormalize", scale.factor = 10000)
hcc29.seu <- FindVariableFeatures(hcc29.seu, selection.method = "vst", nfeatures = 2000)
hcc29_top10 <- head(VariableFeatures(hcc29.seu), 10)
hcc29_all.genes <- rownames(hcc29.seu)
hcc29.seu <- ScaleData(hcc29.seu, features = hcc29_all.genes)

hcc29.seu <- RunPCA(hcc29.seu, verbose = FALSE)
hcc29.seu <- RunUMAP(hcc29.seu, dims = 1:30, verbose = FALSE)
hcc29.seu <- FindNeighbors(hcc29.seu, dims = 1:30, verbose = FALSE)
hcc29.seu <- FindClusters(hcc29.seu, verbose = FALSE)

DimPlot(hcc29.seu, label = TRUE) + NoLegend()
DimPlot(hcc29.seu, label = TRUE,group.by = "orig.ident")
FeaturePlot(hcc29.seu,features = "ITGAM")
FeaturePlot(hcc29.seu,features = "ITGAE")
FeaturePlot(hcc29.seu,features = "ITGAX")
FeaturePlot(hcc29.seu,features =c("AMBP","LMOD1","KRT19","CDH5","CD3D","CD79A","LYZ","COL1A2","ENG"))
VlnPlot(hcc29.seu,features = "GADD45A",group.by = "orig.ident")
DimPlot(hcc28.seu, label = TRUE,group.by = "orig.ident")
DimPlot(hcc29.seu, label = TRUE,group.by = "sample_pt")


hcc28.seu$sample <- "hcc28"
hcc28.seu$sample_pt <- paste(hcc28.seu$sample,hcc28.seu$orig.ident,sep = "_")
hcc29.seu$sample <- "hcc29"
hcc29.seu$sample_pt <- paste(hcc29.seu$sample,hcc29.seu$orig.ident,sep = "_")

hcc.big <- merge(hcc28.seu,y=hcc29.seu)
hcc.big <- PercentageFeatureSet(hcc.big, pattern = "^MT-", col.name = "percent.mt")
hcc.big  <- SCTransform(hcc.big , vars.to.regress = "percent.mt")
hcc.big <- RunPCA(hcc.big, verbose = FALSE)
hcc.big <- RunUMAP(hcc.big, dims = 1:30, verbose = FALSE)
hcc.big  <- FindNeighbors(hcc.big, dims = 1:30, verbose = FALSE)
hcc.big  <- FindClusters(hcc.big, verbose = FALSE)
DimPlot(hcc.big, label = F,group.by = "sample_pt")
FeaturePlot(hcc.big,features =c("AMBP","LMOD1","KRT19","CDH5","CD3D","CD79A","LYZ","COL1A2","ENG"))
DimPlot(hcc.big, label = TRUE)
FeaturePlot(hcc.big,features =c("CD4"))

hcc.big_clusters <- hcc.big$seurat_clusters

hcc.big_celltype = case_when(
  hcc.big_clusters  %in% c("3","19","13","25","30")~"HPC",
  hcc.big_clusters %in% c("29","20","23","6")~"TEC",
  hcc.big_clusters %in% c("17","26")~"B",
  hcc.big_clusters %in% c("21")~"CAF",
  hcc.big_clusters %in% c("22")~"Treg",
  hcc.big_clusters %in% c("1")~"CD4+T",
  hcc.big_clusters %in% c("7","11","2","10","7","14","18")~"CD8+T",
  hcc.big_clusters %in% c("0","5","8","15","28","27")~"My",
  hcc.big_clusters %in% c("9")~"Neutrophil",
  hcc.big_clusters %in% c("16")~"Monocyte",
  hcc.big_clusters %in% c("24")~"DC",
  hcc.big_clusters %in% c("4","12")~"NK",
  TRUE ~ as.character(hcc.big_clusters))
hcc.big@meta.data$celltype= hcc.big_celltype
DimPlot(hcc.big, label = T,group.by = "celltype")
DimPlot(hcc.big, label = T,group.by = "sample")
DimPlot(hcc.big, label = T)


hcc.big <- RunHarmony(hcc.big, "sample_pt")
names(hcc.big@reductions)
hcc.big <- RunUMAP(hcc.big,  dims = 1:15, 
                   reduction = "harmony")
DimPlot(hcc.big,reduction = "umap",label=T,group.by = "sample" ) 
DimPlot(hcc.big,reduction = "umap",label=F,group.by = "sample_pt" ) 
DimPlot(hcc.big,reduction = "umap", label = T,group.by = "celltype")
hcc.big  <- FindClusters(hcc.big, verbose = FALSE)
FeaturePlot(hcc.big,features =c("CD4","CD8A"))
all.markers <- FindAllMarkers(hcc.big)
FeaturePlot(hcc.big,features = "FCN1")
FeaturePlot(hcc.big,features =c("CD4","CD8A"))



cell.prop<-as.data.frame(prop.table(table(hcc.big$sample_pt,hcc.big$celltype)))
colnames(cell.prop) <- c("HCC","origin","proportion")
ggplot(cell.prop,aes(HCC,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 10,color="black"))




immune.big <- subset(hcc.big, subset = celltype %in% c("myeloid","CD4+T","CD8+T","Treg","Neutrophil","Monocyte","DC","NK","B"))


immune_cell.prop<-as.data.frame(prop.table(table(immune.big$sample_pt,immune.big$celltype)))
colnames(immune_cell.prop) <- c("HCC","origin","proportion")
ggplot(immune_cell.prop,aes(HCC,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 10,color="black"))



HPC <- subset(hcc.big, subset = celltype %in% c("HPC"))
HPC <- RunPCA(HPC, verbose = FALSE)
HPC <- RunUMAP(HPC, dims = 1:30, verbose = FALSE)
HPC <- FindNeighbors(HPC, dims = 1:30, verbose = FALSE)
HPC <- FindClusters(HPC, verbose = FALSE)


DimPlot(HPC,group.by = "sample")
HPC <- RunHarmony(HPC, "sample_pt")
HPC <- RunUMAP(HPC,  dims = 1:15, 
                   reduction = "harmony")
DimPlot(HPC, label = F,group.by = "sample_pt")
DimPlot(HPC, label = F,group.by = "sample")

hcc3$orig.ident <- hcc3$sample
hcc3$sample <- "hcc3"
hcc.big_celltype = case_when(
  hcc.big_clusters  %in% c("3","19","13","25","30")~"HPC",
  hcc.big_clusters %in% c("29","20","23","6")~"TEC",
  hcc.big_clusters %in% c("17","26")~"B",
  hcc.big_clusters %in% c("21")~"CAF",
  hcc.big_clusters %in% c("22")~"Treg",
  hcc.big_clusters %in% c("1")~"CD4+T",
  hcc.big_clusters %in% c("7","11","2","10","7","14","18")~"CD8+T",
  hcc.big_clusters %in% c("0","5","8","15","28","27")~"My",
  hcc.big_clusters %in% c("9")~"Neutrophil",
  hcc.big_clusters %in% c("16")~"Monocyte",
  hcc.big_clusters %in% c("24")~"DC",
  hcc.big_clusters %in% c("4","12")~"NK",
  TRUE ~ as.character(hcc.big_clusters))
hcc3$sample_pt <- paste(hcc3$sample,hcc3$orig.ident,sep = "_")
table(hcc3$sample_pt)
hcc_merge <- merge(hcc.big,y=hcc3)

DimPlot(hcc_merge,group.by = "sample")



hcc_merge <- PercentageFeatureSet(hcc_merge, pattern = "^MT-", col.name = "percent.mt")
#hcc_merge  <- SCTransform(hcc_merge , vars.to.regress = "percent.mt")
hcc_merge<- NormalizeData(hcc_merge, normalization.method = "LogNormalize", scale.factor = 10000)
hcc_merge <- FindVariableFeatures(hcc_merge, selection.method = "vst", nfeatures = 2000)
hcc_merge_all.genes <- rownames(hcc_merge)
hcc_merge <- ScaleData(hcc_merge, features = hcc_merge_all.genes)

hcc_merge <- RunPCA(hcc_merge, verbose = FALSE)

hcc_merge  <- FindNeighbors(hcc_merge, dims = 1:30, verbose = FALSE)
hcc_merge  <- FindClusters(hcc_merge,reduction.type = "pca", dims.use = 1:5, resolution = 0.6, print.output = 0, save.SNN = TRUE)
DimPlot(hcc_merge, label = F,group.by = "sample_pt")
DimPlot(hcc_merge, label = F,group.by = "sample")
hcc_merge <- RunUMAP(hcc_merge, dims = 1:30, verbose = FALSE)
hcc_merge <- RunHarmony(hcc_merge, "sample_pt")
hcc_merge <- RunUMAP(hcc_merge,  dims = 1:15, 
                   reduction = "harmony")
DimPlot(hcc_merge,reduction = "umap",label=T) 
DimPlot(hcc_merge,reduction = "umap",label=F,group.by = "sample_pt" ) 
DimPlot(hcc_merge,reduction = "umap", label = T,group.by = "celltype",repel = T)
DimPlot(hcc_merge,reduction = "umap", label = T,group.by = "cell_type")

FeaturePlot(hcc_merge,features = "FOXP3")
FeaturePlot(hcc_merge,features = "FCN1")
FeaturePlot(hcc_merge,features = "S100A9")
FeaturePlot(hcc_merge,features = "S100A8")

hcc2x_celltype <- as.data.frame(hcc_merge$celltype)
hcc3_celltype <- as.data.frame(hcc_merge$cell_type)
hccm_celltype <- merge(hcc2x_celltype,hcc3_celltype,by="row.names")
colnames(hccm_celltype) <- c("row.names","hcc2x","hcc3")
hccm_celltype$hcc2x <- gsub("CD4+T ","CD4+ T Cell",hccm_celltype$hcc2x)  


hcc_merge$celltype[which(hcc_merge$cell_type == "CD4+ T Cell")]<-"CD4+T"
hcc_merge$celltype[which(hcc_merge$cell_type == "CD8+ T Cell")]<-"CD8+T"  
hcc_merge$celltype[which(hcc_merge$cell_type == "B Cell")]<-"B"
hcc_merge$celltype[which(hcc_merge$cell_type == "Treg")]<-"Treg"
hcc_merge$celltype[which(hcc_merge$cell_type == "mki67+ T Cell")]<-"mki67+ T Cell"
hcc_merge$celltype[which(hcc_merge$cell_type == "NK")]<-"NK"
hcc_merge$celltype[which(hcc_merge$cell_type == "myeloid")]<-"Mp"
hcc_merge$celltype[which(hcc_merge$cell_type == "CAF")]<-"CAF"
hcc_merge$celltype[which(hcc_merge$cell_type == "TEC")]<-"TEC"
hcc_merge$celltype[which(hcc_merge$cell_type == "HPC")]<-"HPC"
hcc_merge$celltype[which(hcc_merge$cell_type == "TH2")]<-"CD8+T"
hcc_merge$celltype[which(hcc_merge$celltype == "My")]<-"Mp"
hcc_merge$celltype[which(hcc_merge$cell_type == "Cholangiocytes")]<-"Cholangiocytes"



hcc_merge_clusters <- hcc_merge$seurat_clusters
DimPlot(hcc_merge,label = T,repel = T)
VlnPlot(hcc_merge,features = c("CD4","CD8A"))
FeaturePlot(hcc_merge,features = c("CD4","CD8A"))
DimPlot(hcc_merge,reduction = "umap", label = T,group.by = "celltype",repel = T)

hcc_merge_rough_celltype = case_when(
  hcc_merge_clusters  %in% c("34","4","3","20","23","19","33")~"HPC",
  hcc_merge_clusters %in% c("1","24","36","21")~"TEC",
  hcc_merge_clusters %in% c("16")~"B",
  hcc_merge_clusters %in% c("22")~"CAF",
  hcc_merge_clusters %in% c("10")~"Treg",
  hcc_merge_clusters %in% c("0","14")~"CD4+T",
  hcc_merge_clusters %in% c("12","25","8","5","9","31","18","29")~"CD8+T",
  hcc_merge_clusters %in% c("2","35","7","11")~"My",
  hcc_merge_clusters %in% c("27")~"Neutrophil",
  hcc_merge_clusters %in% c("15")~"Monocyte",
  hcc_merge_clusters %in% c("24")~"DC",
  hcc_merge_clusters %in% c("6","13")~"NK",
  TRUE ~ as.character(hcc_merge_clusters))

hcc_merge@meta.data$rough_celltype= hcc_merge_rough_celltype


DimPlot(hcc_merge,reduction = "umap", label = T,group.by = "rough_celltype",repel = T)




`%notin%` <- Negate(`%in%`)

hcc_merge_anno <- subset(hcc_merge,subset = rough_celltype %notin% c("26","30","28","32","17"))

DimPlot(hcc_merge_anno,reduction = "umap", label = T,group.by = "rough_celltype",repel = T)

DimPlot(hcc_merge_anno,reduction = "umap", label = T,repel = T)


hcc_merge_sample_pt <- hcc_merge$sample_pt

hcc.big_methtype = case_when(
  hcc_merge_sample_pt %in% c("hcc28_NT","hcc29_NT","hcc3_nt")~"normal",
  hcc_merge_sample_pt %in% c("hcc28_PT4","hcc29_PT4","hcc3_pt3","hcc3_pt4")~"hyper",
  hcc_merge_sample_pt %in% c("hcc28_PT1","hcc28_PT2","hcc29_PT1","hcc29_PT3","hcc3_pt1","hcc3_pt2")~"hypo",
  TRUE ~ as.character(hcc_merge_sample_pt))
hcc_merge@meta.data$methtype = hcc.big_methtype
DimPlot(hcc_merge,group.by = "methtype")




hcc_merge_hpc <- subset(hcc_merge,subset=hcc_merge$celltype=="HPC")
DimPlot(hcc_merge_hpc)

hcc_2829_hpc <- subset(hcc.big, subset = celltype == "HPC")
hcc_2829_hpc  <- SCTransform(hcc_2829_hpc , vars.to.regress = "percent.mt")


hcc_2829_hpc <- RunPCA(hcc_2829_hpc, verbose = FALSE)

hcc_2829_hpc  <- FindNeighbors(hcc_2829_hpc, dims = 1:30, verbose = FALSE)
hcc_2829_hpc <- FindClusters(hcc_2829_hpc,reduction.type = "pca", dims.use = 1:5, resolution = 0.6, print.output = 0, save.SNN = TRUE)
DimPlot(hcc_2829_hpc, label = T,group.by = "sample_pt",repel = T)
DimPlot(hcc_merge_hpc, label = F,group.by = "sample")
DimPlot(hcc_merge_hpc, label = F,group.by = "sample")

hcc3_har <- RunHarmony(hcc3, "orig.ident")
hcc3_har <- RunUMAP(hcc3_har,  dims = 1:15, 
                     reduction = "harmony")
DimPlot(hcc3_har,group.by = "orig.ident")
DimPlot(hcc3_har,group.by = "cell_type")


DotPlot(hcc_2829_hpc,features = "GADD45A",group.by = "sample_pt")
DotPlot(hcc_2829_hpc,features = "SNHG6",group.by = "sample_pt")





immune_merge <- subset(hcc_merge_anno, subset = rough_celltype %in% c("B","My","CD4+T","CD8+T","Treg","Neutrophil","Monocyte","NK"))


immune_cell.prop<-as.data.frame(prop.table(table(immune_merge$sample_pt,immune_merge$rough_celltype)))
colnames(immune_cell.prop) <- c("HCC","origin","proportion")
ggplot(immune_cell.prop,aes(HCC,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 10,color="black"))


FeaturePlot(hcc_merge_anno,features = c("AMBP","S100A9","CDH5","CD3D","CD79A","LYZ","COL1A2","ENG","FCN1","CD4","CD8A","FOXP3"))







cell.prop<-as.data.frame(prop.table(table(hcc_merge_anno$sample_pt,hcc_merge_anno$rough_celltype)))
colnames(cell.prop) <- c("HCC","origin","proportion")
ggplot(cell.prop,aes(HCC,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 10,color="black"))


hcc_merge_anno
table(hcc_merge_anno$sample)










hpc_markers_hcc28_4vs1  <- FindMarkers(hcc_2829_hpc,ident.1 = "hcc28_PT4",ident.2 = "hcc28_PT1",group.by = "sample_pt")

hpc_markers_hcc28_4vs1_filt <- subset(hpc_markers_hcc28_4vs1,subset = p_val_adj < 0.05)
hpc_markers_hcc28_4vs1_up <- subset(hpc_markers_hcc28_4vs1_filt,subset = avg_log2FC > 0)
hpc_markers_hcc28_4vs1_down <- subset(hpc_markers_hcc28_4vs1_filt,subset = avg_log2FC < 0)

hpc_markers_hcc28_4vs2  <- FindMarkers(hcc_2829_hpc,ident.1 = "hcc28_PT4",ident.2 = "hcc28_PT2",group.by = "sample_pt")
hpc_markers_hcc28_4vs2_filt <- subset(hpc_markers_hcc28_4vs2,subset = p_val_adj < 0.05)
hpc_markers_hcc28_4vs2_up <- subset(hpc_markers_hcc28_4vs2_filt,subset = avg_log2FC > 0)
hpc_markers_hcc28_4vs2_down <- subset(hpc_markers_hcc28_4vs2_filt,subset = avg_log2FC < 0)

hpc_markers_hcc28_4_common_down <- intersect(row.names(hpc_markers_hcc28_4vs2_down),row.names(hpc_markers_hcc28_4vs1_down))


hpc_markers_hcc28_4vs1_down_gene <- row.names(hpc_markers_hcc28_4vs1_down)
hpc_markers_hcc28_4vs2_down_gene <- row.names(hpc_markers_hcc28_4vs2_down)
hpc_markers_hcc28_4vs1_up_gene <- row.names(hpc_markers_hcc28_4vs1_up)
hpc_markers_hcc28_4vs2_up_gene <- row.names(hpc_markers_hcc28_4vs2_up)
hpc_markers_hcc28_4_common_up <- intersect(hpc_markers_hcc28_4vs1_up_gene,hpc_markers_hcc28_4vs2_up_gene) 


hpc_markers_hcc28_3vs1  <- FindMarkers(hcc_2829_hpc,ident.1 = "hcc29_PT3",ident.2 = "hcc29_PT1",group.by = "sample_pt")
hpc_markers_hcc28_4vs1_filt <- subset(hpc_markers_hcc28_1vs4,subset = p_val_adj < 0.05)
hpc_markers_hcc28_4vs1_up <- subset(hpc_markers_hcc28_1vs4_filt,subset = avg_log2FC > 0)
hpc_markers_hcc28_4vs1_down <- subset(hpc_markers_hcc28_1vs4_filt,subset = avg_log2FC < 0)




hpc_markers_hcc28_4_common_up_go<- enrichGO(gene  = hpc_markers_hcc28_4_common_up,
                                  OrgDb      = org.Hs.eg.db,
                                  keyType    = 'SYMBOL',
                                  ont        = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
hpc_markers_hcc28_4_common_up_go <- as.data.frame(hpc_markers_hcc28_4_common_up_go@result)
hpc_markers_hcc28_4_common_up_go [,"logp"] <- -log10(hpc_markers_hcc28_4_common_up_go $pvalue)
ggplot(data = hpc_markers_hcc28_4_common_up_go[1:20,])+
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
       title="Up Pathways in pt1")



hpc_markers_hcc29_1vs3  <- FindMarkers(hcc_2829_hpc,ident.1 = "hcc29_PT1",ident.2 = "hcc29_PT3",group.by = "sample_pt")
hpc_markers_hcc29_1vs3_filt <- subset(hpc_markers_hcc29_1vs3,subset = p_val_adj < 0.05)
hpc_markers_hcc29_1vs3_up <- subset(hpc_markers_hcc29_1vs3_filt,subset = avg_log2FC > 0)
hpc_markers_hcc29_1vs3_down <- subset(hpc_markers_hcc29_1vs3_filt,subset = avg_log2FC < 0)












hpc_markers_hcc29_1vs3_up_go <- enrichGO(gene  = row.names(hpc_markers_hcc29_1vs3_up),
                                            OrgDb      = org.Hs.eg.db,
                                            keyType    = 'SYMBOL',
                                            ont        = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05)
hpc_markers_hcc29_1vs3_up_go <- as.data.frame(hpc_markers_hcc29_1vs3_up_go@result)
hpc_markers_hcc29_1vs3_up_go [,"logp"] <- -log10(hpc_markers_hcc29_1vs3_up_go $pvalue)
ggplot(data = hpc_markers_hcc29_1vs3_up_go[1:20,])+
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
       title="Up Pathways in pt1")










hpc_markers_hcc29_1vs3_down_go <- enrichGO(gene  = row.names(hpc_markers_hcc29_1vs3_down),
                                         OrgDb      = org.Hs.eg.db,
                                         keyType    = 'SYMBOL',
                                         ont        = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.05)
hpc_markers_hcc29_1vs3_down_go <- as.data.frame(hpc_markers_hcc29_1vs3_down_go@result)
hpc_markers_hcc29_1vs3_down_go [,"logp"] <- -log10(hpc_markers_hcc29_1vs3_down_go $pvalue)
ggplot(data = hpc_markers_hcc29_1vs3_down_go[1:20,])+
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
       title="Down Pathways in pt1")











hcc_merge_anno_T <- subset(hcc_merge,subset= )
