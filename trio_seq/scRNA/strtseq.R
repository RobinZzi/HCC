setwd("~/projects/hcc/analysis/trio_seq/scRNA/strt_sum")
save.image("strt_sum.Rdata")

hcc2_strt <- readRDS("hcc2_strt.RDS")
DimPlot(hcc2_strt,group.by = "sample_pt")
hcc2_strt$sample_pt <- paste("HCC2",hcc2_strt$orig.ident,sep = "_")

hcc3_strt <- readRDS("hcc3_strt.RDS")
hcc3_strt$sample_pt <- paste("HCC3",hcc3_strt$orig.ident,sep = "_")
DimPlot(hcc3_strt,group.by = "sample_pt")

hcc4_strt <- readRDS("hcc4_strt.rds")

hcc4_strt <- PercentageFeatureSet(hcc4_strt, pattern = "^MT-", col.name = "percent.mt")
hcc4_strt <- SCTransform(hcc4_strt, vars.to.regress = "percent.mt", verbose = FALSE)
hcc4_strt <- RunPCA(hcc4_strt, verbose = FALSE)
hcc4_strt <- RunUMAP(hcc4_strt,  dims = 1:15)
hcc4_strt <- FindNeighbors(hcc4_strt, dims = 1:30, verbose = FALSE)
hcc4_strt <- FindClusters(hcc4_strt, verbose = FALSE)
DimPlot(hcc4_strt,group.by = "orig.ident")
hcc4_strt$sample_pt <- paste("HCC4",tolower(hcc4_strt$orig.ident),sep = "_")
DimPlot(hcc4_strt,group.by = "sample_pt")


hcc5_strt <- readRDS("hcc5_strt.rds")

hcc5_strt <- PercentageFeatureSet(hcc5_strt, pattern = "^MT-", col.name = "percent.mt")
hcc5_strt <- SCTransform(hcc5_strt, vars.to.regress = "percent.mt", verbose = FALSE)
hcc5_strt <- RunPCA(hcc5_strt, verbose = FALSE)
hcc5_strt <- RunUMAP(hcc5_strt,  dims = 1:15)
hcc5_strt <- FindNeighbors(hcc5_strt, dims = 1:30, verbose = FALSE)
hcc5_strt <- FindClusters(hcc5_strt, verbose = FALSE)
DimPlot(hcc5_strt,group.by = "orig.ident")
hcc5_strt$sample_pt <- paste("HCC5",tolower(hcc5_strt$orig.ident),sep = "_")
DimPlot(hcc5_strt,group.by = "sample_pt")

hcc6_strt <- readRDS("hcc6_strt.RDS")
DimPlot(hcc6_strt,group.by = "orig.ident")
hcc6_strt$sample_pt <- paste("HCC6",hcc6_strt$orig.ident,sep = "_")
DimPlot(hcc6_strt,group.by = "sample_pt")

hcc7_strt <- readRDS("hcc7_strt.RDS")
DimPlot(hcc7_strt,group.by = "orig.ident")
hcc7_strt$sample_pt <- paste("HCC7",hcc7_strt$orig.ident,sep = "_")
DimPlot(hcc7_strt,group.by = "sample_pt")

hcc8_strt <- readRDS("hcc8_strt.RDS")
DimPlot(hcc8_strt,group.by = "orig.ident")
hcc8_strt$sample_pt <- paste("HCC8",tolower(hcc8_strt$sample_origin),sep = "_")
DimPlot(hcc8_strt,group.by = "sample_pt")

table(hcc8_strt$orig.ident)

hcc9_strt <- readRDS("hcc9_strt.RDS")
DimPlot(hcc9_strt,group.by = "orig.ident")
hcc9_strt$sample_pt <- paste("HCC9",tolower(hcc9_strt$sample_origin),sep = "_")
DimPlot(hcc9_strt,group.by = "sample_pt")

strt_merge <- merge(hcc2_strt,y=c(hcc3_strt,hcc4_strt,hcc5_strt,
                                  hcc6_strt,hcc7_strt,hcc8.seu,hcc9.seu),
                    add.cell.id=c("HCC2","HCC3","HCC4","HCC5","HCC6","HCC7","HCC8","HCC9"))
strt_merge$patient <- str_split_i(strt_merge$sample_pt,"_",1)

strt_merge <- PercentageFeatureSet(strt_merge, pattern = "^MT-", col.name = "percent.mt")
strt_merge <- SCTransform(strt_merge, vars.to.regress = "percent.mt", verbose = FALSE)
strt_merge <- RunPCA(strt_merge, verbose = FALSE)
strt_merge <- RunHarmony(strt_merge, "patient")
strt_merge <- RunUMAP(strt_merge,  dims = 1:15,reduction = "harmony")
strt_merge <- FindNeighbors(strt_merge, dims = 1:10, verbose = FALSE)
strt_merge <- FindClusters(strt_merge, verbose = FALSE)

DimPlot(strt_merge,group.by = "patient")

DimPlot(strt_merge,label = T)

FeaturePlot(strt_merge,features = 'AMBP')   #HPC
FeaturePlot(strt_merge,features = 'COL1A2') #Fibroblast
FeaturePlot(strt_merge,features = 'ENG')    #Endothelial
FeaturePlot(strt_merge,features = 'PTPRC')  #Immune

DimPlot(strt_merge,cells.highlight = "HCC2_pt4_R28")
DimPlot(strt_merge,cells.highlight = "HCC2_pt4_R24")
DimPlot(strt_merge,cells.highlight = "HCC2_pt4_R8")
DimPlot(strt_merge,cells.highlight = "HCC2_pt4_R7")


strt_merge_hpc <- subset(strt_merge,subset=seurat_clusters %in% c(1,2,4,5,6,8,9,10))

hcc2_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC2")
hcc3_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC3")
hcc4_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC4")
hcc5_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC5")
hcc6_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC6")
hcc7_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC7")
hcc8_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC8")
hcc9_strt_hpc <- subset(strt_merge_hpc,subset=patient=="HCC9")

saveRDS(hcc2_strt_hpc,"hcc2_strt_hpc.RDS")
saveRDS(hcc3_strt_hpc,"hcc3_strt_hpc.RDS")
saveRDS(hcc4_strt_hpc,"hcc4_strt_hpc.RDS")
saveRDS(hcc5_strt_hpc,"hcc5_strt_hpc.RDS")
saveRDS(hcc6_strt_hpc,"hcc6_strt_hpc.RDS")
saveRDS(hcc7_strt_hpc,"hcc7_strt_hpc.RDS")
saveRDS(hcc8_strt_hpc,"hcc8_strt_hpc.RDS")
saveRDS(hcc9_strt_hpc,"hcc9_strt_hpc.RDS")
saveRDS(strt_merge,"strt_merge.RDS")



 


hcc8.seu <- CreateSeuratObject(counts = hcc28.strt.raw)
hcc8_ori <- hcc8.seu$orig.ident
hcc8_sample_origin = case_when(
  hcc8_ori %in% c("nt.nt")~"NT",
  hcc8_ori %in% c("pt1a.pt1","pt1b.pt1","pt1c.pt1")~"PT1",
  hcc8_ori %in% c("pt2a.pt2","pt2b.pt2","pt2c.pt2")~"PT2",
  hcc8_ori %in% c("pt4a.pt4","pt4b.pt4","pt4c.pt4")~"PT4",
  TRUE ~ as.character(hcc8_ori))
hcc8.seu@meta.data$sample_origin= hcc8_sample_origin

hcc8.seu$paitent <- 'HCC8'
hcc8.seu$sample_pt <- paste('HCC8',hcc8.seu$sample_origin,sep = "_")


hcc8.seu <- PercentageFeatureSet(hcc8.seu, pattern = "^MT-", col.name = "percent.mt")
hcc8.seu <- SCTransform(hcc8.seu, vars.to.regress = "percent.mt", verbose = FALSE)
hcc8.seu <- RunPCA(hcc8.seu, verbose = FALSE)
hcc8.seu <- RunUMAP(hcc8.seu,  dims = 1:15)
hcc8.seu <- FindNeighbors(hcc8.seu, dims = 1:30, verbose = FALSE)
hcc8.seu <- FindClusters(hcc8.seu, verbose = FALSE)
DimPlot(hcc8.seu,group.by = "orig.ident")
DimPlot(hcc8.seu,group.by = "sample_pt")
saveRDS(hcc8.seu,"hcc8_strt.RDS")


hcc9.seu <- CreateSeuratObject(counts = hcc29.strt.raw)
hcc9.seu$paitent <- 'HCC9'
hcc9_ori <- hcc9.seu$orig.ident
hcc9_sample_origin = case_when(
  hcc9_ori %in% c("nt.nt")~"NT",
  hcc9_ori %in% c("pt1a.pt1","pt1b.pt1","pt1c.pt1")~"PT1",
  hcc9_ori %in% c("pt3a.pt3","pt3b.pt3","pt3c.pt3")~"PT3",
  hcc9_ori %in% c("pt4a.pt4","pt4b.pt4","pt4c.pt4")~"PT4",
  TRUE ~ as.character(hcc9_ori))
hcc9.seu@meta.data$sample_origin= hcc9_sample_origin
hcc9.seu$sample_pt <- paste('HCC9',hcc9.seu$sample_origin,sep = "_")


hcc9.seu <- PercentageFeatureSet(hcc9.seu, pattern = "^MT-", col.name = "percent.mt")
hcc9.seu <- SCTransform(hcc9.seu, vars.to.regress = "percent.mt", verbose = FALSE)
hcc9.seu <- RunPCA(hcc9.seu, verbose = FALSE)
hcc9.seu <- RunUMAP(hcc9.seu,  dims = 1:15)
hcc9.seu <- FindNeighbors(hcc9.seu, dims = 1:30, verbose = FALSE)
hcc9.seu <- FindClusters(hcc9.seu, verbose = FALSE)
DimPlot(hcc9.seu,group.by = "orig.ident")
DimPlot(hcc9.seu,group.by = "sample_pt")
saveRDS(hcc9.seu,"hcc9_strt.RDS")
