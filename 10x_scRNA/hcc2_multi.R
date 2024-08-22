rm(list=ls())
library(harmony)
setwd("~/projects/hcc/analysis/merged_scrna/hcc2_multilib")
HCC2_strt <- readRDS("~/projects/hcc/analysis/merged_scrna/hcc2_multilib/HCC2_strt.RDS")
HCC2_dropseq_PT1 <- readRDS("~/projects/hcc/analysis/merged_scrna/hcc2_multilib/HCC2_dropseq_PT1.rds")
HCC2_10x <- readRDS("~/projects/hcc/analysis/merged_scrna/hcc2_multilib/HCC2_10x.RDS")

save.image("hcc2_merge.Rdata")

HCC2_strt$lib.method <- "strt"
HCC2_dropseq_PT1

hcc2 <- merge(HCC2_strt,y = c(HCC2_dropseq_PT1,HCC2_10x))



hcc2 <- PercentageFeatureSet(hcc2, pattern = "^MT-", col.name = "percent.mt")
hcc2 <- SCTransform(hcc2, vars.to.regress = "percent.mt", verbose = FALSE)
hcc2 <- RunPCA(hcc2, verbose = FALSE)
hcc2 <- RunHarmony(hcc2, "lib.method")
hcc2 <- RunUMAP(hcc2,  dims = 1:15, 
                reduction = "harmony")
hcc2 <- FindNeighbors(hcc2, dims = 1:30, verbose = FALSE)
hcc2 <- FindClusters(hcc2, verbose = FALSE)

DimPlot(hcc2, label = F,group.by = "lib.method",cols= c("#e9e4d4","#57C3F3","#E95C59"),shuffle = T)
DimPlot(hcc2, label = F,group.by = "lib.method",cols= my36colors,shuffle = T)
DimPlot(hcc2, label = T,group.by = "new.ident")+ NoLegend()
DimPlot(hcc2, label = T)+ NoLegend()

hcc2_cluster <- hcc2$seurat_clusters

hcc2_celltype = case_when(
  hcc2_cluster  %in% c("2","24","16","18")~"Endothelial Cell",
  hcc2_cluster %in% c("8", "3","21", "17","14", "22")~"Hepatic Progenitor Cell",
  hcc2_cluster %in% c("11")~"Mast Cell",
  hcc2_cluster %in% c("15")~"Fibroblast Cell",
  hcc2_cluster %in% c("1", "4","26", "5","0", "25","10","6","20")~"T Cell",
  hcc2_cluster %in% c("9")~"NK Cell",
  hcc2_cluster %in% c("27", "23","7", "19","12")~"Myeloid Cell",
  hcc2_cluster %in% c("28", "13")~"B Cell",
  TRUE ~ as.character(hcc2_cluster))

hcc2$cell_type <- hcc2_celltype

DimPlot(hcc2, label = T,group.by = "cell_type",cols=my36colors)+ NoLegend()


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175','#FFEEAD')#颜色设置



VlnPlot(hcc2,features = c("nFeature_RNA", "nCount_RNA"),group.by = "lib.method",pt.size = 0,cols= c("#e9e4d4","#57C3F3","#E95C59"))




hcc2_hpc <- subset(hcc2,subset=cell_type=="Hepatic Progenitor Cell")

hcc2_hpc <- PercentageFeatureSet(hcc2_hpc, pattern = "^MT-", col.name = "percent.mt")
hcc2_hpc <- SCTransform(hcc2_hpc, vars.to.regress = "percent.mt", verbose = FALSE)
hcc2_hpc <- RunPCA(hcc2_hpc, verbose = FALSE)
hcc2_hpc <- RunHarmony(hcc2_hpc, "lib.method")
hcc2_hpc <- RunUMAP(hcc2_hpc,  dims = 1:15,reduction = "harmony")
hcc2_hpc <- FindNeighbors(hcc2_hpc, dims = 1:30, verbose = FALSE)
hcc2_hpc <- FindClusters(hcc2_hpc, verbose = FALSE)

DimPlot(hcc2_hpc, label = T,group.by = "patient_pt",cols=my36colors)+ NoLegend()
DimPlot(hcc2_hpc, label = T,group.by = "lib.method",cols=my36colors)+ NoLegend()
DimPlot(hcc2_hpc, label = T,group.by = "patient",cols=my36colors)+ NoLegend()
DimPlot(hcc2_hpc, label = T,group.by = "orig.ident",cols=my36colors)+ NoLegend()


hcc2_hpc_PT1 <- subset(hcc2_hpc , subset = orig.ident == "PT1")
hcc2_PT1_distance <- as.data.frame(Embeddings(object = hcc2_hpc_PT1[["umap"]]))
hcc2_hpc_PT2 <- subset(hcc2_hpc , subset = orig.ident == "PT2")
hcc2_PT2_distance <- as.data.frame(Embeddings(object = hcc2_hpc_PT2[["umap"]]))
hcc2_hpc_PT3 <- subset(hcc2_hpc , subset = orig.ident == "PT3")
hcc2_PT3_distance <- as.data.frame(Embeddings(object = hcc2_hpc_PT3[["umap"]]))


hcc3_distance <- as.data.frame(rbind(colMeans(hcc2_PT1_distance),colMeans(hcc2_PT2_distance),colMeans(hcc2_PT3_distance)))
row.names(hcc3_distance) <- c("PT1","PT2","PT3")


dist_hcc3 = dist(hcc3_distance, method = "euclidean")
hclust_dist_hcc3 = hclust(dist_hcc3, method = "complete")
plot(hclust_dist_hcc3)



hcc2_hpc_sub <- subset(hcc2_hpc_sub,subset = orig.ident %in% c('PT1','PT2','PT3'))
hcc2_hpc_sub <- PercentageFeatureSet(hcc2_hpc_sub, pattern = "^MT-", col.name = "percent.mt")
hcc2_hpc_sub <- SCTransform(hcc2_hpc_sub, vars.to.regress = "percent.mt", verbose = FALSE)
hcc2_hpc_sub <- RunPCA(hcc2_hpc_sub, verbose = FALSE)
hcc2_hpc_sub <- RunHarmony(hcc2_hpc_sub, "lib.method")
hcc2_hpc_sub <- RunUMAP(hcc2_hpc_sub,  dims = 1:15,reduction = "harmony")
hcc2_hpc_sub <- FindNeighbors(hcc2_hpc_sub, dims = 1:30, verbose = FALSE)
hcc2_hpc_sub <- FindClusters(hcc2_hpc_sub, verbose = FALSE)
DimPlot(hcc2_hpc_sub, label = F,group.by = "orig.ident",cols=c('#8C549C','#D6E7A3','#E59CC4'))

        