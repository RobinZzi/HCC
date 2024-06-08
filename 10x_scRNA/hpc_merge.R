rm(list = ls())
library(stringr)
library(Seurat)
library(dplyr)
library(harmony)
library(clusterProfiler)
library(DoubletFinder)
library(glmGamPoi)
library(ggtree)
library(phylogram)
setwd("~/projects/hcc/data/10x_scRNA/hcc3-10x/result")
setwd("~/projects/hcc/analysis/10x_scRNA/")
saveRDS(hcc_big_HPC,"hpc.RDS")

save.image("hpc.Rdata")
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

hcc28_cnt <- readRDS("~/projects/hcc/data/10x_scRNA/hcc28.10x.raw.cnt.rds")
hcc28_meta <- readRDS("~/projects/hcc/data/10x_scRNA/hcc28.10x.meta.rds")
hcc29_cnt <- readRDS("~/projects/hcc/data/10x_scRNA/hcc29.10x.raw.cnt.rds")
hcc29_meta <- readRDS("~/projects/hcc/data/10x_scRNA/hcc29.10x.meta.rds")


hcc3_meta$sample <- "hcc3"
hcc28_meta$sample <- "hcc28"
hcc29_meta$sample <- "hcc29"

hcc3_meta$sample_pt <- paste("hcc3",hcc3_meta$orig.ident,sep = "_")  
hcc28_meta$sample_pt <- paste("hcc28",hcc28_meta$orig.ident,sep = "_")
hcc29_meta$sample_pt <- paste("hcc29",hcc29_meta$orig.ident,sep = "_")  

hcc3_sce <- CreateSeuratObject(counts = hcc3_cnt,meta.data = hcc3_meta)

hcc28_sce <- CreateSeuratObject(counts = hcc28_cnt,meta.data = hcc28_meta)

hcc29_sce <- CreateSeuratObject(counts = hcc29_cnt,meta.data = hcc29_meta)

VlnPlot(hcc3_sce,features = "GADD45A",group.by = "sample_pt")
VlnPlot(hcc28_sce,features = "TNNC1",group.by = "sample_pt")

hcc_big <- merge(hcc3_sce,y=hcc28_sce)
hcc_big <- merge(hcc_big,y=hcc29_sce)


hcc_big[["percent.mt"]] <- PercentageFeatureSet(hcc_big, pattern = "^MT-")

HB.genes <- c("HBA1","HBA2")
HB_m <- match(HB.genes, rownames(hcc_big@assays$RNA)) 
HB.genes <- rownames(hcc_big@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
hcc_big[["percent.HB"]]<-PercentageFeatureSet(hcc_big, features=HB.genes) 


col.num <- length(levels(as.factor(hcc_big@meta.data$sample_pt)))



VlnPlot(hcc_big, group.by = "sample_pt",  
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0.00, 
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())




hcc_big_filt <- subset(hcc_big, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA < 20000)

VlnPlot(hcc_big_filt, group.by = "sample_pt",  
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0.00, 
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())



hcc_big_filt <- SCTransform(hcc_big_filt,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
hcc_big_filt<- RunPCA(hcc_big_filt, verbose = FALSE)
hcc_big_filt <- RunUMAP(hcc_big_filt, dims = 1:30, verbose = FALSE)

hcc_big_filt <- FindNeighbors(hcc_big_filt, dims = 1:30, verbose = FALSE)
hcc_big_filt <- FindClusters(hcc_big_filt, verbose = FALSE)
DimPlot(hcc_big_filt, label = TRUE) + NoLegend()
DimPlot(hcc_big_filt, label = F,group.by = "sample")
hcc_big_filt <- RunHarmony(hcc_big_filt, "sample_pt")
hcc_big_filt <- RunUMAP(hcc_big_filt,  dims = 1:15, reduction = "harmony")
DimPlot(hcc_big_filt,group.by = "sample")
DimPlot(hcc_big_filt,label = T)

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

DimPlot(hcc29_strt,group.by = "sample_origin")


hcc_big_filt_removedb <- hcc_big_filt


hcc_big_filt_removedb <- subset(hcc_big_filt,subset = DF.classifications_0.25_0.19_4281 == "Singlet")



DimPlot(hcc_big_filt_removedb,group.by = "sample")
DimPlot(hcc_big_filt_removedb,label = T)
hcc_big_filt_removedb <- subset(hcc_big_filt_removedb,subset= seurat_clusters != 34)

FeaturePlot(hcc_big_filt_removedb,features = "FCGR3A")
FeaturePlot(hcc_big_filt_removedb,features = "MKI67")
FeaturePlot(hcc_big_filt,features = c("AMBP","S100A9","CDH5","CD3D","CD79A","LYZ","LMOD1","ENG","FCN1","CD4","CD8A","FOXP3"))
hcc_merge_clusters <- hcc_big_filt_removedb$seurat_clusters


hcc_merge_rough_celltype = case_when(
  hcc_merge_clusters  %in% c("31","35","8","6","15","34")~"HPC",
  TRUE ~ as.character(hcc_merge_clusters))

hcc_big_filt_removedb@meta.data$rough_celltype= hcc_merge_rough_celltype

hcc_big_filt_removedb$sample_pt

DimPlot(hcc_big_filt_removedb,group.by = "rough_celltype",label = T)
DimPlot(hcc_big_filt_removedb,group.by = "sample",label = T)
DimPlot(hcc_big_filt_removedb,label = T)

VlnPlot(hcc_big_filt_removedb,features = "FCGR3A")

hcc_big_HPC <- subset(hcc_big_filt_removedb,subset = rough_celltype == "HPC")


hcc_big_HPC_big <- merge(hcc_big_HPC,y=hcc29strt)

hcc_big_HPC <- RunPCA(hcc_big_HPC, verbose = FALSE)
hcc_big_HPC <- RunUMAP(hcc_big_HPC, dims = 1:30, verbose = FALSE)

hcc_big_HPC <- FindNeighbors(hcc_big_HPC, dims = 1:30, verbose = FALSE)
hcc_big_HPC <- FindClusters(hcc_big_HPC, verbose = FALSE)
DimPlot(hcc_big_HPC,group.by = "tech",label = T,repel = T)


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










hcc28_strt <- readRDS("~/projects/hcc/data/10x_scRNA/hcc3-10x/result/hcc28-strt.RDS")
hcc28_strt$sample_pt <- paste("hcc28",hcc28_strt$sample_origin,sep = "_")
hcc29_strt <- readRDS("~/projects/hcc/data/10x_scRNA/hcc3-10x/result/hcc29-strt.RDS")
hcc29_strt$sample_pt <- paste("hcc29",hcc29_strt$sample_origin,sep = "_")


hcc_big_HPC$tech <- "10x"
hcc28_strt$tech <- "strt"
hcc29_strt$tech <- "strt"
hcc28_strt$sample <- "hcc28"
hcc29_strt$sample <- "hcc29"


DimPlot(hcc28_strt,group.by = "sample_origin")


saveRDS(hcc29_strt,"hcc29_strt.RDS")

saveRDS(hcc28_strt,"hcc28_strt.RDS")








hcc3_strt_cnt <- readRDS("~/projects/hcc/data/10x_scRNA/result/raw.count.rds")
fs <- list.files("./",pattern = "*.txt")
for (i in 1:length(fs)) {
  data <- fread(fs[i],col.names = c("barcode","sample"))
  if(i == 1){
    hcc3_strt_meta <- data
  }else{
    hcc3_strt_meta <- rbind(hcc3_strt_meta,data)
  }
}

hcc3_strt_meta <- as.data.frame(hcc3_strt_meta)
hcc3_strt_meta[,3:4] <- as.data.frame(str_split_fixed(hcc3_strt_meta$sample,"_",2))
row.names(hcc3_strt_meta) <- hcc3_strt_meta$barcode



hcc3_strt_meta <- as.data.frame(colnames(hcc3_strt_cnt))
hcc3_strt_meta[,2:3] <- as.data.frame(str_split_fixed(hcc3_strt_meta[,1],"t.",2))
colnames(hcc3_strt_meta) <- c("a","b","c")
hcc3_strt_meta <- select(hcc3_strt_meta,a,c)
hcc3_strt_meta[,3:4] <- as.data.frame(str_split_fixed(hcc3_strt_meta$c,"_",2))
row.names(hcc3_strt_meta) <- hcc3_strt_meta$a
hcc3_strt_meta <- select(hcc3_strt_meta,V1)
colnames(hcc3_strt_meta) <- "sample"
hcc3_strt_meta$sample <- gsub(".pt1","pt1",hcc3_strt_meta$sample)
hcc3_strt_meta$sample <- gsub(".pt2","pt2",hcc3_strt_meta$sample)
hcc3_strt_meta$sample <- gsub(".pt3","pt3",hcc3_strt_meta$sample)
hcc3_strt_meta$sample <- gsub(".pt4","pt4",hcc3_strt_meta$sample)
hcc3_strt_meta$sample <- toupper(hcc3_strt_meta$sample)
hcc3_strt_meta$sample_pt <- paste("hcc3",hcc3_strt_meta$sample,sep = "_")
hcc3_strt_meta$tech <- "strt"
hcc3_strt_meta$sample <- "hcc3"
hcc3_strt <- CreateSeuratObject(counts = hcc3_strt_cnt,meta.data = hcc3_strt_meta)
hcc3_strt[["percent.mt"]] <- PercentageFeatureSet(hcc3_strt, pattern = "^MT-")


VlnPlot(hcc3_strt,features = c("percent.mt","nCount_RNA","nFeature_RNA"),group.by = "sample_pt")

hpc_merge <- merge(hcc_big_HPC,y=hcc28_strt)
hpc_merge <- merge(hpc_merge,y=hcc29_strt)

hpc_merge <- merge(hpc_merge,y=hcc3_strt)
hpc_merge <- SCTransform(hpc_merge,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)



hpc_merge <- RunPCA(hpc_merge, verbose = FALSE)
hpc_merge <- RunUMAP(hpc_merge, dims = 1:30, verbose = FALSE)

hpc_merge <- FindNeighbors(hpc_merge, dims = 1:30, verbose = FALSE)
hpc_merge <- FindClusters(hpc_merge, verbose = FALSE)
DimPlot(hpc_merge,group.by = "sample_pt",label = T,repel = T)
DimPlot(hpc_merge,group.by = "tech",label = T,repel = T)
hpc_merge_harmony <- RunHarmony(hpc_merge, "tech")
hpc_merge_harmony  <- RunUMAP(hpc_merge_harmony ,  dims = 1:15, reduction = "harmony")
DimPlot(hpc_merge_harmony,group.by = "tech",label = T,repel = T)


DimPlot(hpc_merge_harmony,group.by = "sample_pt",split.by = "tech",label = T,repel = T)
DimPlot(hpc_merge_harmony,group.by = "sample",split.by = "tech",label = T,repel = T)
DimPlot(hpc_merge_harmony,group.by = "sample_pt",label = T,repel = T)

FeaturePlot(hpc_merge_harmony,features = c("hcc_PT1","hcc3_pt1"))



 

hpc_sample_origin <- hpc_merge_harmony$sample_pt



hpc_meth_type = case_when(
  hpc_sample_origin  %in% c("hcc3_PT3","hcc3_PT4")~"hcc3_admeth",
  hpc_sample_origin  %in% c("hcc3_PT1","hcc3_PT2")~"hcc3_demeth",
  hpc_sample_origin  %in% c("hcc28_PT1","hcc28_PT2","hcc28_PT4")~"hcc28_demeth",
  hpc_sample_origin  %in% c("hcc29_PT1","hcc29_PT3","hcc29_PT4")~"hcc29_demeth",
  hpc_sample_origin  %in% c("hcc3_NT","hcc28_NT","hcc29_NT")~"NT",
  TRUE ~ as.character(hpc_sample_origin))


hpc_merge_harmony$hpc_meth_type <- hpc_meth_type







hcc3_demeth <- FindMarkers(hpc_merge_harmony,ident.1 = "hcc3_demeth",ident.2 = "hcc3_admeth",group.by = "hpc_meth_type")
hcc28_demeth <- FindMarkers(hpc_merge_harmony,ident.1 = "hcc28_demeth",ident.2 = "hcc3_admeth",group.by = "hpc_meth_type")
hcc29_demeth <- FindMarkers(hpc_merge_harmony,ident.1 = "hcc29_demeth",ident.2 = "hcc3_admeth",group.by = "hpc_meth_type")



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

up_venn <- list(hcc3_demeth_up$gene,hcc28_demeth_up$gene,hcc29_demeth_up$gene)


venn.diagram(up_venn, filename = 'up.png', imagetype = 'png', ,category.names = c("hcc3" , "hcc28" , "hcc29"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

down_venn <- list(hcc3_demeth_down$gene,hcc28_demeth_down$gene,hcc29_demeth_down$gene)

venn.diagram(down_venn, filename = 'down.png', imagetype = 'png', ,category.names = c("hcc3" , "hcc28" , "hcc29"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')


write.table(hcc3_demeth_up,"hcc3_demeth_up.txt")
write.table(hcc3_demeth_down,"hcc3_demeth_down.txt")
write.table(hcc28_demeth_up,"hcc28_demeth_up.txt")
write.table(hcc28_demeth_down,"hcc28_demeth_down.txt")
write.table(hcc29_demeth_up,"hcc29_demeth_up.txt")
write.table(hcc29_demeth_down,"hcc29_demeth_down.txt")






VlnPlot(hcc_big_HPC,features = "GADD45A",group.by = "sample_pt")
VlnPlot(hcc_big_HPC,features = "SNHG6",group.by = "sample_pt")



hcc3_hpc_infer <- subset(hcc_big_HPC, subset= sample =="hcc3")
hcc28_hpc_infer <- subset(hcc_big_HPC, subset= sample =="hcc28")
hcc29_hpc_infer <- subset(hcc_big_HPC, subset= sample =="hcc29")
DimPlot(hcc3_hpc_infer,group.by = "sample_pt",label=T)
DimPlot(hcc28_hpc_infer,group.by = "sample_pt",label=T,repel = T)+NoLegend()
DimPlot(hcc29_hpc_infer,group.by = "sample_pt",label=T)

saveRDS(hcc3_hpc_infer,"hcc3_hpc.RDS")
saveRDS(hcc28_hpc_infer,"hcc28_hpc.RDS")
saveRDS(hcc29_hpc_infer,"hcc29_hpc.RDS")


hcc3_hpc_infer <- RunPCA(hcc3_hpc_infer, verbose = FALSE)
hcc3_hpc_infer <- RunUMAP(hcc3_hpc_infer, dims = 1:30, verbose = FALSE)
hcc3_hpc_infer <- FindNeighbors(hcc3_hpc_infer, dims = 1:30, verbose = FALSE)
hcc3_hpc_infer <- FindClusters(hcc3_hpc_infer, verbose = FALSE)
DimPlot(hcc3_hpc_infer,group.by = "sample_pt",label = T,label.size = 6)+NoLegend()



hcc3_NT <- subset(hcc3_hpc_infer , subset = sample_pt == "hcc3_NT")
hcc3_NT_distance <- as.data.frame(Embeddings(object = hcc3_NT[["umap"]]))
hcc3_PT1 <- subset(hcc3_hpc_infer , subset = sample_pt == "hcc3_PT1")
hcc3_PT1_distance <- as.data.frame(Embeddings(object = hcc3_PT1[["umap"]]))
hcc3_PT2 <- subset(hcc3_hpc_infer , subset = sample_pt == "hcc3_PT2")
hcc3_PT2_distance <- as.data.frame(Embeddings(object = hcc3_PT2[["umap"]]))
hcc3_PT3 <- subset(hcc3_hpc_infer , subset = sample_pt == "hcc3_PT3")
hcc3_PT3_distance <- as.data.frame(Embeddings(object = hcc3_PT3[["umap"]]))
hcc3_PT4 <- subset(hcc3_hpc_infer , subset = sample_pt == "hcc3_PT4")
hcc3_PT4_distance <- as.data.frame(Embeddings(object = hcc3_PT4[["umap"]]))
hcc3_distance <- as.data.frame(rbind(colMeans(hcc3_NT_distance),colMeans(hcc3_PT1_distance),colMeans(hcc3_PT2_distance),
                                     colMeans(hcc3_PT3_distance),colMeans(hcc3_PT4_distance)))
row.names(hcc3_distance) <- c("NT","PT1","PT2","PT3","PT4")

hcc3_distance <- hcc3_distance[-1,]

dist_hcc3 = dist(hcc3_distance, method = "euclidean")
hclust_dist_hcc3 = hclust(dist_hcc3, method = "complete")
plot(hclust_dist_hcc3)
plot(hclust_dist_hcc28)



hcc28_hpc_infer <- RunPCA(hcc28_hpc_infer, verbose = FALSE)
hcc28_hpc_infer <- RunUMAP(hcc28_hpc_infer, dims = 1:30, verbose = FALSE)
hcc28_hpc_infer <- FindNeighbors(hcc28_hpc_infer, dims = 1:30, verbose = FALSE)
hcc28_hpc_infer <- FindClusters(hcc28_hpc_infer, verbose = FALSE)
DimPlot(hcc28_hpc_infer,group.by = "sample_pt")


hcc28_NT <- subset(hcc28_hpc_infer , subset = sample_pt == "hcc28_NT")
hcc28_NT_distance <- as.data.frame(Embeddings(object = hcc28_NT[["umap"]]))
hcc28_PT1 <- subset(hcc28_hpc_infer , subset = sample_pt == "hcc28_PT1")
hcc28_PT1_distance <- as.data.frame(Embeddings(object = hcc28_PT1[["umap"]]))
hcc28_PT2 <- subset(hcc28_hpc_infer , subset = sample_pt == "hcc28_PT2")
hcc28_PT2_distance <- as.data.frame(Embeddings(object = hcc28_PT2[["umap"]]))
hcc28_PT4 <- subset(hcc28_hpc_infer , subset = sample_pt == "hcc28_PT4")
hcc28_PT4_distance <- as.data.frame(Embeddings(object = hcc28_PT4[["umap"]]))
hcc28_distance <- as.data.frame(rbind(colMeans(hcc28_NT_distance),colMeans(hcc28_PT1_distance),
                                      colMeans(hcc28_PT2_distance),colMeans(hcc28_PT4_distance)))
row.names(hcc28_distance) <- c("NT","PT1","PT2","PT4")

hcc28_distance <- hcc28_distance[-1,]
dist_hcc28 = dist(hcc28_distance, method = "euclidean")
hclust_dist_hcc28 = hclust(dist_hcc28, method = "complete")
plot(hclust_dist_hcc28)




hcc29_hpc_infer <- RunPCA(hcc29_hpc_infer, verbose = FALSE)
hcc29_hpc_infer <- RunUMAP(hcc29_hpc_infer, dims = 1:30, verbose = FALSE)
hcc29_hpc_infer <- FindNeighbors(hcc29_hpc_infer, dims = 1:30, verbose = FALSE)
hcc29_hpc_infer <- FindClusters(hcc29_hpc_infer, verbose = FALSE)
DimPlot(hcc29_hpc_infer,group.by = "sample_pt")



hcc29_NT <- subset(hcc29_hpc_infer , subset = sample_pt == "hcc29_NT")
hcc29_NT_distance <- as.data.frame(Embeddings(object = hcc29_NT[["umap"]]))
hcc29_PT1 <- subset(hcc29_hpc_infer , subset = sample_pt == "hcc29_PT1")
hcc29_PT1_distance <- as.data.frame(Embeddings(object = hcc29_PT1[["umap"]]))
hcc29_PT3 <- subset(hcc29_hpc_infer , subset = sample_pt == "hcc29_PT3")
hcc29_PT3_distance <- as.data.frame(Embeddings(object = hcc29_PT3[["umap"]]))
hcc29_PT4 <- subset(hcc29_hpc_infer , subset = sample_pt == "hcc29_PT4")
hcc29_PT4_distance <- as.data.frame(Embeddings(object = hcc29_PT4[["umap"]]))
hcc29_distance <- as.data.frame(rbind(colMeans(hcc29_NT_distance),colMeans(hcc29_PT1_distance),
                                      colMeans(hcc29_PT3_distance),colMeans(hcc29_PT4_distance)))
row.names(hcc29_distance) <- c("NT","PT1","PT3","PT4")

hcc29_distance <- hcc29_distance[-1,]
dist_hcc29 = dist(hcc29_distance, method = "euclidean")
hclust_dist_hcc29 = hclust(dist_hcc29, method = "complete")
plot(hclust_dist_hcc29)



hcc3_PT1_pca <- as.data.frame(Embeddings(hcc3_PT1, reduction = "pca"))
hcc3_PT2_pca <- as.data.frame(Embeddings(hcc3_PT2, reduction = "pca"))
hcc3_PT3_pca <- as.data.frame(Embeddings(hcc3_PT3, reduction = "pca"))
hcc3_PT4_pca <- as.data.frame(Embeddings(hcc3_PT4, reduction = "pca"))
hcc3_pca <- as.data.frame(rbind(colMeans(hcc3_PT1_pca),colMeans(hcc3_PT2_pca),colMeans(hcc3_PT3_pca),colMeans(hcc3_PT4_pca)))
row.names(hcc3_pca) <- c("PT1","PT2","PT3","PT4")
pca_hcc3 = dist(hcc3_pca, method = "euclidean")
hclust_pca_hcc3 = hclust(pca_hcc3, method = "complete")
plot(hclust_pca_hcc3)
plot(hclust_pca_hcc28)
hclust_pca_hcc3_tree <-cutree(hclust_pca_hcc3,3)
hcc3_d = data.frame(label=names(hclust_pca_hcc3_tree),member=factor(hclust_pca_hcc3_tree))
ggtree(as.phylo(hclust_pca_hcc3), linetype='dashed', color = "#487AA1") %<+% hcc3_d +
  geom_tiplab(aes(color = member),size=4) +
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+NoLegend()
  

ggtree(as.phylo(hclust_dist_hcc29),layout="rectangular",branch.length="none", color = "#487AA1")+
  geom_tiplab(hjust =0.1,size=6,fontface="plain")+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

ggtree(as.phylo(hclust_pca_hcc28),layout="rectangular",branch.length="none", color = "#487AA1")+
  geom_tiplab(hjust =0.1,size=6,fontface="plain")+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

ggtree(as.phylo(hclust_pca_hcc29),layout="rectangular",branch.length="none", color = "#487AA1")+
  geom_tiplab(hjust =0.1,size=6,fontface="plain")+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

hcc28_PT1_pca <- as.data.frame(Embeddings(hcc28_PT1, reduction = "pca"))
hcc28_PT2_pca <- as.data.frame(Embeddings(hcc28_PT2, reduction = "pca"))
hcc28_PT4_pca <- as.data.frame(Embeddings(hcc28_PT4, reduction = "pca"))
hcc28_pca <- as.data.frame(rbind(colMeans(hcc28_PT1_pca),colMeans(hcc28_PT2_pca),colMeans(hcc28_PT4_pca)))
row.names(hcc28_pca) <- c("PT1","PT2","PT4")
pca_hcc28 = dist(hcc28_pca, method = "euclidean")
hclust_pca_hcc28 = hclust(pca_hcc28, method = "complete")
plot(hclust_pca_hcc28)






hcc29_PT1_pca <- as.data.frame(Embeddings(hcc29_PT1, reduction = "pca"))
hcc29_PT3_pca <- as.data.frame(Embeddings(hcc29_PT3, reduction = "pca"))
hcc29_PT4_pca <- as.data.frame(Embeddings(hcc29_PT4, reduction = "pca"))
hcc29_pca <- as.data.frame(rbind(colMeans(hcc29_PT1_pca),colMeans(hcc29_PT3_pca),colMeans(hcc29_PT4_pca)))
row.names(hcc29_pca) <- c("PT1","PT3","PT4")
pca_hcc29 = dist(hcc29_pca, method = "euclidean")
hclust_pca_hcc29 = hclust(pca_hcc29, method = "complete")
plot(hclust_pca_hcc29)
  


hcc_big_HPC_sample_pt <- hcc_big_HPC$sample_pt


hcc_big_HPC_methtype = case_when(
  hcc_big_HPC_sample_pt  %in% c("hcc28_NT","hcc29_NT","hcc3_NT","hcc3_PT3","hcc3_PT4")~"Consistently-methylated Samples",
  hcc_big_HPC_sample_pt %in% c("hcc28_PT1","hcc28_PT2","hcc28_PT4","hcc29_PT1","hcc29_PT4","hcc29_PT3","hcc3_PT1","hcc3_PT2")~"De-methylated Samples",
  TRUE ~ as.character(hcc_big_HPC_sample_pt))
hcc_big_HPC$methtype <- hcc_big_HPC_methtype
DimPlot(hcc_big_HPC,group.by = "methtype",label=F)+NoLegend()
DimPlot(hcc_big_HPC,group.by = "sample",label=T,label.size = 8)+NoLegend()
VlnPlot(hcc_big_HPC,features = "SNHG6",group.by = "sample",split.by = "methtype")
VlnPlot(hcc_big_HPC,features = "TNNC1",group.by = "sample",split.by = "methtype")
VlnPlot(hcc_big_HPC,features = "GADD45A",group.by = "sample",split.by = "methtype")
