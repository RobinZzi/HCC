rm(list = ls())
library(stringr)
library(Seurat)
library(dplyr)
library(harmony)
library(clusterProfiler)
library(DoubletFinder)
setwd("~/projects/hcc/data/10x_scRNA/hcc3-10x/result")

hcc3_cnt <- Read10X("~/projects/hcc/data/10x_scRNA/hcc3-10x/result/hcc3")
hcc3_meta <- fread("barcode_sample.txt")
hcc3_meta_info <- as.data.frame(colnames(hcc3_cnt))
hcc3_meta_info[,2:3] <- str_split_fixed(hcc3_meta_info$`colnames(hcc3_cnt)`,"-",2) 
colnames(hcc3_meta_info) <- c("barcode","cell_barcode","sample_id")

hcc3_meta_info$sample_id[hcc3_meta_info$sample_id== 1] <- "NT"
hcc3_meta_info$sample_id[hcc3_meta_info$sample_id== 2] <- "PT1"
hcc3_meta_info$sample_id[hcc3_meta_info$sample_id== 3] <- "PT2"
hcc3_meta_info$sample_id[hcc3_meta_info$sample_id== 4] <- "PT3"
hcc3_meta_info$sample_id[hcc3_meta_info$sample_id== 5] <- "PT4"

hcc3_meta <- select(hcc3_meta_info,sample_id)
row.names(hcc3_meta) <- hcc3_meta_info$barcode
colnames(hcc3_meta) <- "orig.ident"

hcc3_meta$sample <- "hcc3"
hcc28_meta$sample <- "hcc28"
hcc29_meta$sample <- "hcc29"

hcc3_meta$sample_pt <- paste("hcc3",hcc3_meta$orig.ident,sep = "_")  
hcc28_meta$sample_pt <- paste("hcc28",hcc28_meta$orig.ident,sep = "_")
hcc29_meta$sample_pt <- paste("hcc29",hcc29_meta$orig.ident,sep = "_")  

hcc3_sce <- CreateSeuratObject(counts = hcc3_cnt,meta.data = hcc3_meta)

hcc28_sce <- CreateSeuratObject(counts = hcc28_cnt,meta.data = hcc28_meta)

hcc29_sce <- CreateSeuratObject(counts = hcc29_cnt,meta.data = hcc29_meta)




hcc_big <- merge(hcc3_sce,y=hcc28_sce)
hcc_big <- merge(hcc_big,y=hcc29_sce)


hcc_big[["percent.mt"]] <- PercentageFeatureSet(hcc_big, pattern = "^MT-")

HB.genes <- c("HBA1","HBA2")
HB_m <- match(HB.genes, rownames(hcc_big@assays$RNA)) 
HB.genes <- rownames(hcc_big@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
hcc_big[["percent.HB"]]<-PercentageFeatureSet(hcc_big, features=HB.genes) 


col.num <- length(levels(as.factor(hcc_big@meta.data$sample)))



VlnPlot(hcc_big_filt, group.by = "sample",  
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0.01, 
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())




hcc_big_filt <- subset(hcc_big, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA < 20000)





hcc_big_filt <- SCTransform(hcc_big_filt, vars.to.regress = "percent.mt", verbose = FALSE)
hcc_big_filt<- RunPCA(hcc_big_filt, verbose = FALSE)
hcc_big_filt <- RunUMAP(hcc_big_filt, dims = 1:30, verbose = FALSE)

hcc_big_filt <- FindNeighbors(hcc_big_filt, dims = 1:30, verbose = FALSE)
hcc_big_filt <- FindClusters(hcc_big_filt, verbose = FALSE)
DimPlot(hcc_big_filt, label = TRUE) + NoLegend()
DimPlot(hcc_big_filt, label = F,group.by = "sample")
hcc_big_filt <- RunHarmony(hcc_big_filt, "sample_pt")
hcc_big_filt <- RunUMAP(hcc_big_filt,  dims = 1:15, reduction = "harmony")


rep_list <- paramSweep_v3(hcc_big_filt,PCs = 1:30,sct=T)
sweep_stats <- summarizeSweep(rep_list,GT=F)
bcmvn <- find.pK(sweep_stats)
opt_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
print(opt_pK)

annotations <- hcc_big_filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.06*nrow(hcc_big_filt@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

hcc_big_filt <- doubletFinder_v3(hcc_big_filt,
                                 PCs = 1:30,
                                 pN = 0.25,
                                 pK = opt_pK,
                                 nExp = nExp_poi,
                                 reuse.pANN = F,
                                 sct = T)


DimPlot(hcc_big_filt,group.by = "DF.classifications_0.25_0.19_4281")







hcc_big_filt_removedb <- subset(hcc_big_filt,subset = DF.classifications_0.25_0.19_4281 == "Singlet")



DimPlot(hcc_big_filt_removedb,group.by = "sample")
DimPlot(hcc_big_filt_removedb,label = T)
hcc_big_filt_removedb <- subset(hcc_big_filt_removedb,subset= seurat_clusters != 34)

FeaturePlot(hcc_big_filt_removedb,features = "FCGR3A")
FeaturePlot(hcc_big_filt_removedb,features = "MKI67")
FeaturePlot(hcc_big_filt_removedb,features = c("AMBP","S100A9","CDH5","CD3D","CD79A","LYZ","LMOD1","ENG","FCN1","CD4","CD8A","FOXP3"))
hcc_merge_clusters <- hcc_big_filt_removedb$seurat_clusters


hcc_merge_rough_celltype = case_when(
  hcc_merge_clusters  %in% c("14","6","22","13","36","32")~"HPC",
  hcc_merge_clusters %in% c("20","3","16","29","27")~"TEC",
  hcc_merge_clusters %in% c("15")~"B",
  hcc_merge_clusters %in% c("17")~"CAF",
  hcc_merge_clusters %in% c("17","21","18","10","30","8","2","5","12")~"T",
  hcc_merge_clusters %in% c("0","1","11","35","28","33","4","31","9")~"My",
  hcc_merge_clusters %in% c("24")~"Neutrophil",
  hcc_merge_clusters %in% c("23","25")~"Monocyte",
  hcc_merge_clusters %in% c("7","19")~"NK",
  hcc_merge_clusters %in% c("26")~"MKI67+ T",
  TRUE ~ as.character(hcc_merge_clusters))

hcc_big_filt_removedb@meta.data$rough_celltype= hcc_merge_rough_celltype



DimPlot(hcc_big_filt_removedb,group.by = "rough_celltype",label = T)
DimPlot(hcc_big_filt_removedb,label = T)

VlnPlot(hcc_big_filt_removedb,features = "FCGR3A")

hcc_big_HPC <- subset(hcc_big_filt_removedb,subset = rough_celltype == "HPC")




hcc_big_HPC <- RunPCA(hcc_big_HPC, verbose = FALSE)
hcc_big_HPC <- RunUMAP(hcc_big_HPC, dims = 1:30, verbose = FALSE)

hcc_big_HPC <- FindNeighbors(hcc_big_HPC, dims = 1:30, verbose = FALSE)
hcc_big_HPC <- FindClusters(hcc_big_HPC, verbose = FALSE)
DimPlot(hcc_big_HPC,group.by = "sample_pt",label = T,repel = T)
hcc_big_HPC <- RunHarmony(hcc_big_HPC, "sample_pt")
hcc_big_HPC<- RunUMAP(hcc_big_HPC,  dims = 1:15, reduction = "harmony")

DimPlot(hcc_big_HPC, group.by = "sample_pt",label = T,repel = T)

hpc_sample_origin <- hcc_big_HPC$sample_pt


hpc_meth_type = case_when(
  hpc_sample_origin  %in% c("hcc3_PT3","hcc3_PT4")~"hcc3_admeth",
  hpc_sample_origin  %in% c("hcc3_PT1","hcc3_PT2")~"hcc3_demeth",
  hpc_sample_origin  %in% c("hcc28_PT1","hcc28_PT2","hcc28_PT4")~"hcc28_demeth",
  hpc_sample_origin  %in% c("hcc29_PT1","hcc29_PT3","hcc29_PT4")~"hcc29_demeth",
  hpc_sample_origin  %in% c("hcc3_NT","hcc28_NT","hcc29_NT")~"NT",
  TRUE ~ as.character(hpc_sample_origin))

hcc_big_HPC$sample_methtype <- hpc_meth_type

DimPlot(hcc_big_HPC, group.by = "sample_methtype",label = T,repel = T)
VlnPlot(hcc_big_HPC,features = "ACTB",group.by = "sample_pt")
VlnPlot(hcc_big_HPC,features = "CTNNB1",group.by = "sample_pt")

DotPlot(hcc_big_HPC,features = "LRP1",group.by = "sample_pt")

hcc3_demeth <- FindMarkers(hcc_big_HPC,ident.1 = "hcc3_demeth",ident.2 = "hcc3_admeth",group.by = "sample_methtype")
hcc28_demeth <- FindMarkers(hcc_big_HPC,ident.1 = "hcc28_demeth",ident.2 = "hcc3_admeth",group.by = "sample_methtype")
hcc29_demeth <- FindMarkers(hcc_big_HPC,ident.1 = "hcc29_demeth",ident.2 = "hcc3_admeth",group.by = "sample_methtype")



hcc3_demeth$gene <- row.names(hcc3_demeth) 
hcc28_demeth$gene <- row.names(hcc28_demeth) 
hcc29_demeth$gene <- row.names(hcc29_demeth) 



hcc3_demeth_up <- subset(hcc3_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc3_demeth_down <- subset(hcc3_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)

hcc28_demeth_up <- subset(hcc28_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc28_demeth_down <- subset(hcc28_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)

hcc29_demeth_up <- subset(hcc29_demeth, subset = p_val_adj < 0.05 & avg_log2FC > 0)
hcc29_demeth_down <- subset(hcc29_demeth, subset = p_val_adj < 0.05 & avg_log2FC < 0)




common_down <- intersect(hcc3_demeth_down$gene,hcc28_demeth_down$gene)
common_down <- intersect(common_down ,hcc29_demeth_down$gene)


common_up <- intersect(hcc3_demeth_up$gene,hcc28_demeth_up$gene)
common_up <- intersect(common_up,hcc29_demeth_up$gene)









hpc_common_up_go <- enrichGO(gene  = common_up,
                                  OrgDb      = org.Hs.eg.db,
                                  keyType    = 'SYMBOL',
                                  ont        = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
hpc_common_up_go <- as.data.frame(hpc_common_up_go@result)
hpc_common_up_go [,"logp"] <- -log10(hpc_common_up_go$pvalue)
ggplot(data = hpc_common_up_go[1:20,])+
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








hpc_common_down_go <- enrichGO(gene  = common_down,
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
hpc_common_down_go <- as.data.frame(hpc_common_down_go@result)
hpc_common_down_go [,"logp"] <- -log10(hpc_common_down_go$pvalue)
ggplot(data = hpc_common_down_go[1:20,])+
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

VlnPlot(hcc_big_filt_removedb,group.by = "rough_celltype",split.by = "sample",features = "CD4")
