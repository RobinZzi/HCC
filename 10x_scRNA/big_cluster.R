library(Seurat)
library(DoubletFinder)
my37colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175','#FFEEAD')#颜色设置

DimPlot(bigseu,group.by = "new.ident",cols = my36colors)#cell_type
DimPlot(bigseu,cols = my37colors,group.by = "seurat_clusters")#cell_type

cell_cluster_info <- bigseu$new.ident
big_cell_cluster_info = case_when(
  cell_cluster_info  %in% c("Fibroblast")~"Fibroblast",
  cell_cluster_info  %in% c("Endothelial cell")~"Endothelial cell",
  cell_cluster_info  %in% c("HPC")~"Hepatic Progenitor Cell",
  cell_cluster_info  %in% c("Mast cell")~"Mast cell",
  cell_cluster_info  %in% c("PlasmaB cell", "B cell",
                            "Proliferative T","CD8+ memory", "CD8+ exhausted",
                            "CD8+ cytotoxic","CD4+ memory","CD4+ Treg","NK")~"Lymphocyte",
  cell_cluster_info  %in% c("Dendritic cell","Macrophage")~"Myeloid Cell",
  cell_cluster_info  %in% c("Neutrophil")~"Granulocytes",
  TRUE ~ as.character(cell_cluster_info))
bigseu$big_cell_cluster_info <- big_cell_cluster_info
DimPlot(bigseu,group.by = "big_cell_cluster_info",label=F,cols=my37colors)

big_cell_cluster2_info = case_when(
  cell_cluster_info  %in% c("Fibroblast")~"Interstitial Cell",
  cell_cluster_info  %in% c("Endothelial cell")~"Endothelial Cell",
  cell_cluster_info  %in% c("HPC")~"Epithelial Cell",
  cell_cluster_info  %in% c("PlasmaB cell", "B cell",
                            "Proliferative T","CD8+ memory", "CD8+ exhausted",
                            "CD8+ cytotoxic","CD4+ memory","CD4+ Treg","NK",
                            "Dendritic cell","Macrophage","Mast cell",
                            "Neutrophil" )~"Immune Cell",
  TRUE ~ as.character(cell_cluster_info))
bigseu$big_cell_cluster2_info <- big_cell_cluster2_info
DimPlot(bigseu,group.by = "big_cell_cluster2_info",label=F,cols=my37colors[3:6])
FeaturePlot(bigseu,features = "CD68")
FeaturePlot(bigseu,features = "KIT")


FeaturePlot(bigseu,features = c("PTPRC","KRT8","PECAM1","COL1A2"),ncol=2)
