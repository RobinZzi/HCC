my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3',
               '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A',
               '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
               '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35',
               '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD',
               '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')



library(Seurat)
library(dplyr)
save.image("gaqc.Rdata")
# 1. 样本对应的细胞数
cell_counts <- sce.big@meta.data %>%
  group_by(sample) %>%
  summarise(cell_count = n())

# 2. 每个细胞的reads（假设为UMI数，存于nCount_RNA）
# 计算每个样本的平均reads数
mean_reads_per_sample <- sce.big@meta.data %>%
  group_by(sample) %>%
  summarise(mean_reads = mean(nCount_RNA, na.rm = TRUE))

# 3. 总基因数（每个细胞检测到的基因数，存于nFeature_RNA）
total_genes_per_cell <- sce.big@meta.data %>%
  group_by(sample) %>%
  summarise(mean_genes = mean(nFeature_RNA, na.rm = TRUE))

# 4. 中位数基因数
median_genes_per_cell <- sce.big@meta.data %>%
  group_by(sample) %>%
  summarise(median_genes = median(nFeature_RNA, na.rm = TRUE))



# 将所有结果合并
sample_summary <- cell_counts %>%
  left_join(mean_reads_per_sample, by = "sample") %>%
  left_join(total_genes_per_cell, by = "sample") %>%
  left_join(median_genes_per_cell, by = "sample")

print(sample_summary)


sample_summary_tx <- sample_summary
sample_summary_tx$lib <- "10x" 

sample_summary_dp <- sample_summary
sample_summary_dp$lib <- "dropseq" 

sample_summary_st <- sample_summary
sample_summary_st$lib <- "strtseq" 



sce.big <- SCTransform(sce.big,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
sce.big <- RunPCA(sce.big, verbose = FALSE)
sce.big <- RunUMAP(sce.big, dims = 1:30, verbose = FALSE)
DimPlot(sce.big,label = T)
FeaturePlot(sce.big,features = c("AMBP","COL1A2","PTPRC","ENG"))
DimPlot(sce.big,group.by = "patient" )
sce.big <- subset(sce.big,subset=seurat_clusters %in% c(1,4,5,7,8,3,2,0,10,6))
DimPlot(sce.big)
VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "patient")










sce.big <- SCTransform(sce.big,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
sce.big <- RunPCA(sce.big, verbose = FALSE)
sce.big <- RunUMAP(sce.big, dims = 1:30, verbose = FALSE)
DimPlot(sce.big,label = T)
FeaturePlot(sce.big,features = c("AMBP","COL1A2","PTPRC","ENG"))
DimPlot(sce.big,group.by = "patient" )
sce.big <- subset(sce.big,subset=seurat_clusters %in% c(1,4,5,7,8,3,2,0,10,6))
DimPlot(sce.big,label = T)
cell_cluster_info <- sce.big$seurat_clusters
big_cell_cluster_info = case_when(
  cell_cluster_info  %in% c(2)~"Fibroblast",
  cell_cluster_info  %in% c(6,10,0)~"Endothelial cell",
  cell_cluster_info  %in% c(1,8,7,4,5,3)~"Hepatic Progenitor Cell",
  TRUE ~ as.character(cell_cluster_info))
sce.big$cell_type <- big_cell_cluster_info
DimPlot(sce.big,group.by = "cell_type",cols=c(my36colors[3],my36colors[6],my36colors[4]),shuffle = T)

st_patient_info <- sce.big$patient

st_patient_info2 <- case_when(
  st_patient_info  %in% c("hcc11")~"HCC7",
  st_patient_info  %in% c("hcc7")~"HCC6",
  st_patient_info  %in% c("hcc4")~"HCC3",
  st_patient_info  %in% c("hcc3")~"HCC2",
  TRUE ~ as.character(st_patient_info)
)
sce.big$patient2 <- st_patient_info2
VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "patient2",pt.size = 0.00,cols=c('HCC1'='#E95C59','HCC2'='#E5D2DD','HCC3'='#BD956A',
                                                                                                                                'HCC4'='#57C3F3','HCC5'='#3A6963','HCC6'='#F1BB72',
                                                                                                                                'HCC7'='#8C549C','HCC8'='#D6E7A3','HCC9'='#E59CC4'))


bigseu_tx <- SCTransform(bigseu_tx,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
bigseu_tx <- RunPCA(bigseu_tx, verbose = FALSE)
bigseu_tx <- RunUMAP(bigseu_tx, dims = 1:30, verbose = FALSE)
DimPlot(bigseu_tx,label = T)
FeaturePlot(bigseu_tx,features = c("AMBP","COL1A2","PTPRC","ENG"))
DimPlot(bigseu_tx,group.by = "patient" )
DimPlot(bigseu_tx,group.by = "big_cell_cluster2_info",cols = my36colors[3:6])
VlnPlot(bigseu_tx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "patient",pt.size = 0.00,cols=c('HCC1'='#E95C59','HCC2'='#E5D2DD','HCC3'='#BD956A',
                                                                                                                                'HCC4'='#57C3F3','HCC5'='#3A6963','HCC6'='#F1BB72',
                                                                                                                                'HCC7'='#8C549C','HCC8'='#D6E7A3','HCC9'='#E59CC4'))




bigseu_dp <- SCTransform(bigseu_dp,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
bigseu_dp <- RunPCA(bigseu_dp, verbose = FALSE)
bigseu_dp <- RunUMAP(bigseu_dp, dims = 1:30, verbose = FALSE)
DimPlot(bigseu_dp,label = T)
FeaturePlot(bigseu_dp,features = c("AMBP","COL1A2","PTPRC","ENG"))
DimPlot(bigseu_dp,group.by = "patient" )
DimPlot(bigseu_dp,group.by = "big_cell_cluster2_info" )
DimPlot(bigseu_dp,group.by = "new.ident" )
DimPlot(bigseu_dp,group.by = "big_cell_cluster2_info",cols = my36colors[3:6])
VlnPlot(bigseu_dp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "patient",pt.size = 0.00,cols=c('HCC1'='#E95C59','HCC2'='#E5D2DD','HCC3'='#BD956A',
                                                                                                                                 'HCC4'='#57C3F3','HCC5'='#3A6963','HCC6'='#F1BB72',
                                                                                                                                 'HCC7'='#8C549C','HCC8'='#D6E7A3','HCC9'='#E59CC4'))





hpc_tx <- subset(bigseu_tx,subset=big_cell_cluster2_info=="Epithelial Cell")
hpc_dp <- subset(bigseu_dp,subset=big_cell_cluster2_info=="Epithelial Cell")
hpc_st <- subset(sce.big,subset=cell_type=="Hepatic Progenitor Cell")


hpc_merge <- merge(hpc_tx,y=c(hpc_dp,hpc_st))
saveRDS(hpc_merge,"hpc_merge_ga.Rds")


rm(list = ls())
