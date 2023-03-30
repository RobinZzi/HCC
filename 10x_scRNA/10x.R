setwd("~/projects/hcc/single_cell_hcc28_hcc29/10x")
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
rm(list=ls())





save.image("10x.Rdata")

hcc29.seu <- CreateSeuratObject(counts = hcc29_cnt,
                                meta.data = hcc29_meta)
hcc29.seu <- PercentageFeatureSet(hcc29.seu, pattern = "^MT-", col.name = "percent.mt")
hcc29.seu <- subset(hcc29.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
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


table(HPC$sample_pt)
