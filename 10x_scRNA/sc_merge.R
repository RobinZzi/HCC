rm(list = ls())
library(ggplot2)
library(ggpubr)
library(infercnv)
library(devtools)
library(rjags)
library(clusterProfiler)
library(gridExtra)
setwd("~/projects/hcc/analysis/merged_scrna")

save.image("merged_scRNA.Rdata")

devtools::install_github("broadinstitute/infercnv")
bigseu <- readRDS("~/projects/hcc/data/10x_scRNA/merge/HCC-scRNA-seq-inte.rds")

bigseu$lib.method
hcc.big_patients <- bigseu$patient
hcc.big_cntype = case_when(
  hcc.big_patients  %in% c("HCC3","HCC5","HCC7")~"SN",
  hcc.big_patients %in% c("HCC1", "HCC2","HCC4", "HCC6","HCC8", "HCC9")~"CMN",
  TRUE ~ as.character(hcc.big_patients))

bigseu$cntype <- hcc.big_cntype

DimPlot(bigseu,group.by = "group",cols = my36colors[1:3],shuffle = T) #NT,PT,ST

DimPlot(bigseu,group.by = "seurat_clusters",label = T,cols=my37colors)+NoLegend()

DimPlot(bigseu,group.by = "patient",cols=my37colors,shuffle = T)

DimPlot(bigseu,group.by = "lib.method",cols = my36colors[9:10])#10x,dropseq

DimPlot(bigseu,group.by = "orig.ident")# PT1/PT2/PT3

DimPlot(bigseu,group.by = "new.ident",cols = my36colors)#cell_type

DimPlot(bigseu,group.by = "cntype") #CMN,SN

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
bigseu<- CellCycleScoring(bigseu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
table(bigseu$Phase)
DimPlot(bigseu,group.by = "Phase",cols = my37colors,label=F)


VlnPlot(bigseu, features = c("nFeature_RNA", "nCount_RNA"), group.by = "lib.method",cols = my37colors,pt.size = 0)

# 重新调整图层顺序
bigseu@meta.data$patient <- factor(bigseu@meta.data$patient, levels = c( "HCC1", "HCC2", "HCC3",
                                                                      "HCC4", "HCC5", "HCC6",
                                                                      "HCC7","HCC8", "HCC9"))
bigseu@meta.data$new.ident <- factor(bigseu@meta.data$new.ident, levels = rev(c( "HPC", "Fibroblast", "Endothelial cell",
                                                                         "Neutrophil", "PlasmaB cell", "B cell",
                                                                         "Proliferative T","CD8+ memory", "CD8+ exhausted",
                                                                         "CD8+ cytotoxic","CD4+ memory","CD4+ Treg","Mast cell",
                                                                         "NK","Dendritic cell","Macrophage")))
markers <- c('AMBP','COL1A2', 'ENG', 'S100A8', 'S100A9', 'CD79A', 'MS4A1','LTB',
             'MKI67','CD8A','PDCD1', 'CCL5','GZMA', 'CD3D', 'IL32','CD4','FOXP3','SPON2',
              "FCN1","LYZ")

DimPlot(bigseu,group.by = "patient",cols = t(my36colors[1:9]))#patient

HPC <- subset(bigseu,subset=new.ident=="HPC")

FeaturePlot(bigseu,features = "PDCD1")
VlnPlot(bigseu,features = "PDCD1",group.by = "group")
VlnPlot(bigseu,features = "PDCD1",group.by = "new.ident")
DimPlot(HPC,group.by = "patient")

hpc_merge <- SCTransform(HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
hpc_merge <- RunPCA(hpc_merge, verbose = FALSE)
hpc_merge <- RunUMAP(hpc_merge, dims = 1:30, verbose = FALSE)

DimPlot(hpc_merge,group.by = "patient",label=T,repel = T)
DimPlot(hpc_merge,group.by = "group",label=T)
DimPlot(hpc_merge,group.by = "lib.method",label=T)
DimPlot(hpc_merge,group.by = "new.ident",label=F,cols = my36colors )
DimPlot(hpc_merge,group.by = "group",label=T)

cell.prop<-as.data.frame(prop.table(table(bigseu$patient,bigseu$new.ident)))
colnames(cell.prop) <- c("HCC","origin","proportion")
ggplot(cell.prop,aes(HCC,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 10,color="black"))

cell_prop_ptst<-as.data.frame(prop.table(table(bigseu$patient,bigseu$new.ident)))
cell.prop<-as.data.frame(prop.table(table(bigseu$patient,bigseu$new.ident)))

markers <- c('CD3D','CCL5', 'NKG7', 'GZMA', 'IL32', 'CD4','CD8A','FOXP3', 'CD3E', 'LTB', 'S100A8', 'S100A9', 'CD79A', "ENG" ,
             'FCN1', 'MS4A1', 'SPON2','FCER1A','SERPINF1', "LYZ","AMBP","COL1A2","MKI67","PDCD1")
ex.markers <- c('LAYN','PDCD1','CTLA4','HAVCR2')
my37colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175','#FFEEAD')#颜色设置

VlnPlot(bigseu, features = markers,pt.size=0, stack = T,
        cols = my36colors,group.by = "new.ident")+NoLegend()+
        theme(axis.text.y = element_text(face="bold",size = 10,color = 'black'),axis.text.x = element_blank(),
               strip.text.x.top   = element_text(face="bold",size = 10,angle=60,color = 'black'),
              strip.background = element_blank(),strip.clip="off")

DotPlot(bigseu, features = markers,
        group.by = "new.ident")+NoLegend()+
  theme_bw()+theme(axis.text.y = element_text(face="bold",size = 10,color = 'black'),
                   axis.text.x = element_text(face="bold",size = 10,angle=45,hjust=1,color = 'black'))

 VlnPlot(bigseu, features = ex.markers,pt.size=0, stack = T,
        cols = my36colors,group.by = "new.ident")+NoLegend()+
  theme(axis.text.x = element_blank())

VlnPlot(bigseu, features = c("LAYN","HAVCR2"),pt.size=0, stack = T,
        cols = my36colors,group.by = "new.ident")+NoLegend()+
  theme(axis.text.x = element_blank())



VlnPlot(bigseu, features = c("MKI67", "TOP2A", "HAVCR2","LAG3"),pt.size=0, stack = T,
        cols = my36colors,group.by = "new.ident")+NoLegend()+
  theme(axis.text.x = element_blank())
DimPlot(bigseu,group.by = "new.ident",label = T)


HPC_patient_info <- hpc_merge$patient
HPC_subtype_info = case_when(
  HPC_patient_info  %in% c("HCC3","HCC5","HCC7")~"SN",
  HPC_patient_info  %in% c("HCC1","HCC2","HCC4","HCC6","HCC8","HCC9")~"CMN",
  TRUE ~ HPC_patient_info(HPC_patient_info))
hpc_merge$subtype_info <- HPC_subtype_info
DimPlot(hpc_merge,group.by = "subtype_info",label=T)

HCC1_HPC <- subset(HPC,subset=patient=="HCC1")
HCC2_HPC <- subset(HPC,subset=patient=="HCC2")
HCC3_HPC <- subset(HPC,subset=patient=="HCC3")
HCC4_HPC <- subset(HPC,subset=patient=="HCC4")
HCC5_HPC <- subset(HPC,subset=patient=="HCC5")
HCC6_HPC <- subset(HPC,subset=patient=="HCC6")
HCC7_HPC <- subset(HPC,subset=patient=="HCC7")
HCC8_HPC <- subset(HPC,subset=patient=="HCC8")
HCC9_HPC <- subset(HPC,subset=patient=="HCC9")











HCC1_HPC <- SCTransform(HCC1_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC1_HPC <- RunPCA(HCC1_HPC, verbose = FALSE)
HCC1_HPC <- RunUMAP(HCC1_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC1_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC1_HPC,group.by = "group",label=T)
DimPlot(HCC1_HPC,group.by = "satellite_region",label=T)

hcc1_sample_origin <- HCC1_HPC$orig.ident
satellite_region_hcc1 = case_when(
  hcc1_sample_origin  %in% c("PT2","PT4")~"primary",
  hcc1_sample_origin  %in% c("PT5")~"satellite",
  hcc1_sample_origin  %in% c("NT")~"nt",
  TRUE ~ as.character(hcc1_sample_origin))
HCC1_HPC$satellite_region <- satellite_region_hcc1

hcc1_sta_markers <- FindMarkers(HCC1_HPC,ident.1 = "primary",ident.2 = "satellite",group.by = "satellite_region" )
hcc1_sta_markers_sig <- subset(hcc1_sta_markers,subset = p_val_adj < 0.05)
hcc1_sta_markers_sig_up <- subset(hcc1_sta_markers_sig,subset = avg_log2FC > 0)
hcc1_sta_markers_sig_down <- subset(hcc1_sta_markers_sig,subset = avg_log2FC < 0)

hcc1_sta_up_go <- enrichGO(gene  = row.names(hcc1_sta_markers_sig_up),
                         OrgDb      = org.Hs.eg.db,
                         keyType    = 'SYMBOL',
                         ont        = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
hcc1_sta_up_go <- as.data.frame(hcc1_sta_up_go@result)
hcc1_sta_up_go [,"logp"] <- -log10(hcc1_sta_up_go$pvalue)
ggplot(data = hcc1_sta_up_go[1:20,])+
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
       title="hcc1_sta_up")

hcc1_sta_down_go <- enrichGO(gene  = row.names(hcc1_sta_markers_sig_down),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
hcc1_sta_down_go <- as.data.frame(hcc1_sta_down_go@result)
hcc1_sta_down_go [,"logp"] <- -log10(hcc1_sta_down_go$pvalue)
ggplot(data = hcc1_sta_down_go[1:20,])+
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
       title="hcc1_sta_down")

HCC2_HPC <- SCTransform(HCC2_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC2_HPC <- RunPCA(HCC2_HPC, verbose = FALSE)
HCC2_HPC <- RunUMAP(HCC2_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC2_HPC,group.by = "orig.ident",label=T,repel = T)
DimPlot(HCC2_HPC,group.by = "group",label=T)




HCC3_HPC <- SCTransform(HCC3_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC3_HPC <- RunPCA(HCC3_HPC, verbose = FALSE)
HCC3_HPC <- RunUMAP(HCC3_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC3_HPC,group.by = "orig.ident",label=T,repel = T)
DimPlot(HCC3_HPC,group.by = "group",label=T)
hcc3_sample_origin <- HCC3_HPC$orig.ident
satellite_region_hcc3 = case_when(
  hcc3_sample_origin  %in% c("PT1","PT2")~"primary",
  hcc3_sample_origin  %in% c("PT3")~"satellite",
  hcc3_sample_origin  %in% c("NT")~"nt",
  TRUE ~ as.character(hcc3_sample_origin))
HCC3_HPC$satellite_region <- satellite_region_hcc3

hcc3_sta_markers <- FindMarkers(HCC3_HPC,ident.1 = "primary",ident.2 = "satellite",group.by = "satellite_region" )
hcc3_sta_markers_sig <- subset(hcc3_sta_markers,subset = p_val_adj < 0.05)
hcc3_sta_markers_sig_up <- subset(hcc3_sta_markers_sig,subset = avg_log2FC > 0)
hcc3_sta_markers_sig_down <- subset(hcc3_sta_markers_sig,subset = avg_log2FC < 0)

hcc3_sta_up_go <- enrichGO(gene  = row.names(hcc3_sta_markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
hcc3_sta_up_go <- as.data.frame(hcc3_sta_up_go@result)
hcc3_sta_up_go [,"logp"] <- -log10(hcc3_sta_up_go$pvalue)
ggplot(data = hcc3_sta_up_go[1:20,])+
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
       title="hcc3_sta_up")

hcc3_sta_down_go <- enrichGO(gene  = row.names(hcc3_sta_markers_sig_down),
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
hcc3_sta_down_go <- as.data.frame(hcc3_sta_down_go@result)
hcc3_sta_down_go [,"logp"] <- -log10(hcc3_sta_down_go$pvalue)
ggplot(data = hcc3_sta_down_go[1:20,])+
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
       title="hcc3_sta_down")


HCC4_HPC <- SCTransform(HCC4_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC4_HPC <- RunPCA(HCC4_HPC, verbose = FALSE)
HCC4_HPC <- RunUMAP(HCC4_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC4_HPC,group.by = "orig.ident",label=T,repel = T)
DimPlot(HCC4_HPC,group.by = "group",label=T)




HCC5_HPC <- SCTransform(HCC5_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC5_HPC <- RunPCA(HCC5_HPC, verbose = FALSE)
HCC5_HPC <- RunUMAP(HCC5_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC5_HPC,group.by = "orig.ident",label=T,repel = T)
DimPlot(HCC5_HPC,group.by = "group",label=T)

hcc5_sample_origin <- HCC5_HPC$orig.ident
satellite_region_hcc5 = case_when(
  hcc5_sample_origin  %in% c("PT1","PT2","PT3","PT4")~"primary",
  hcc5_sample_origin  %in% c("PT5")~"satellite",
  hcc5_sample_origin  %in% c("NT")~"nt",
  TRUE ~ as.character(hcc5_sample_origin))
HCC5_HPC$satellite_region <- satellite_region_hcc5

hcc5_sta_markers <- FindMarkers(HCC5_HPC,ident.1 = "primary",ident.2 = "satellite",group.by = "satellite_region" )
hcc5_sta_markers_sig <- subset(hcc5_sta_markers,subset = p_val_adj < 0.05)
hcc5_sta_markers_sig_up <- subset(hcc5_sta_markers_sig,subset = avg_log2FC > 0)
hcc5_sta_markers_sig_down <- subset(hcc5_sta_markers_sig,subset = avg_log2FC < 0)

hcc5_sta_up_go <- enrichGO(gene  = row.names(hcc5_sta_markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
hcc5_sta_up_go <- as.data.frame(hcc5_sta_up_go@result)
hcc5_sta_up_go [,"logp"] <- -log10(hcc5_sta_up_go$pvalue)
ggplot(data = hcc5_sta_up_go[1:20,])+
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
       title="hcc5_sta_up")

hcc5_sta_down_go <- enrichGO(gene  = row.names(hcc5_sta_markers_sig_down),
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
hcc5_sta_down_go <- as.data.frame(hcc5_sta_down_go@result)
hcc5_sta_down_go [,"logp"] <- -log10(hcc5_sta_down_go$pvalue)
ggplot(data = hcc5_sta_down_go[1:20,])+
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
       title="hcc5_sta_down")


HCC6_HPC <- SCTransform(HCC6_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC6_HPC <- RunPCA(HCC6_HPC, verbose = FALSE)
HCC6_HPC <- RunUMAP(HCC6_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC6_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC6_HPC,group.by = "group",label=T)

hcc6_sample_origin <- HCC6_HPC$orig.ident
satellite_region_hcc6 = case_when(
  hcc6_sample_origin  %in% c("PT1","PT2","PT3")~"primary",
  hcc6_sample_origin  %in% c("PT4","PT5","PT6")~"satellite",
  hcc6_sample_origin  %in% c("NT")~"nt",
  TRUE ~ as.character(hcc6_sample_origin))
HCC6_HPC$satellite_region <- satellite_region_hcc6

hcc6_sta_markers <- FindMarkers(HCC6_HPC,ident.1 = "primary",ident.2 = "satellite",group.by = "satellite_region" )
hcc6_sta_markers_sig <- subset(hcc6_sta_markers,subset = p_val_adj < 0.05)
hcc6_sta_markers_sig_up <- subset(hcc6_sta_markers_sig,subset = avg_log2FC > 0)
hcc6_sta_markers_sig_down <- subset(hcc6_sta_markers_sig,subset = avg_log2FC < 0)

hcc6_sta_up_go <- enrichGO(gene  = row.names(hcc6_sta_markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
hcc6_sta_up_go <- as.data.frame(hcc6_sta_up_go@result)
hcc6_sta_up_go [,"logp"] <- -log10(hcc6_sta_up_go$pvalue)
ggplot(data = hcc6_sta_up_go[1:20,])+
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
       title="hcc6_sta_up")

hcc6_sta_down_go <- enrichGO(gene  = row.names(hcc6_sta_markers_sig_down),
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
hcc6_sta_down_go <- as.data.frame(hcc6_sta_down_go@result)
hcc6_sta_down_go [,"logp"] <- -log10(hcc6_sta_down_go$pvalue)
ggplot(data = hcc6_sta_down_go[1:20,])+
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
       title="hcc6_sta_down")


HCC7_HPC <- SCTransform(HCC7_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC7_HPC <- RunPCA(HCC7_HPC, verbose = FALSE)
HCC7_HPC <- RunUMAP(HCC7_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC7_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC7_HPC,group.by = "group",label=T)


hcc7_sample_origin <- HCC7_HPC$orig.ident
satellite_region_hcc7 = case_when(
  hcc7_sample_origin  %in% c("PT1","PT2","PT3","PT4","PT5")~"primary",
  hcc7_sample_origin  %in% c("PT6")~"satellite",
  hcc7_sample_origin  %in% c("NT")~"nt",
  TRUE ~ as.character(hcc7_sample_origin))
HCC7_HPC$satellite_region <- satellite_region_hcc7

hcc7_sta_markers <- FindMarkers(HCC7_HPC,ident.1 = "primary",ident.2 = "satellite",group.by = "satellite_region" )
hcc7_sta_markers_sig <- subset(hcc7_sta_markers,subset = p_val_adj < 0.05)
hcc7_sta_markers_sig_up <- subset(hcc7_sta_markers_sig,subset = avg_log2FC > 0)
hcc7_sta_markers_sig_down <- subset(hcc7_sta_markers_sig,subset = avg_log2FC < 0)

hcc7_sta_up_go <- enrichGO(gene  = row.names(hcc7_sta_markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
hcc7_sta_up_go <- as.data.frame(hcc7_sta_up_go@result)
hcc7_sta_up_go [,"logp"] <- -log10(hcc7_sta_up_go$pvalue)
ggplot(data = hcc7_sta_up_go[1:20,])+
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
       title="hcc7_sta_up")

hcc7_sta_down_go <- enrichGO(gene  = row.names(hcc7_sta_markers_sig_down),
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
hcc7_sta_down_go <- as.data.frame(hcc7_sta_down_go@result)
hcc7_sta_down_go [,"logp"] <- -log10(hcc7_sta_down_go$pvalue)
ggplot(data = hcc7_sta_down_go[1:20,])+
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
       title="hcc7_sta_down")


HCC8_HPC <- SCTransform(HCC8_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC8_HPC <- RunPCA(HCC8_HPC, verbose = FALSE)
HCC8_HPC <- RunUMAP(HCC8_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC8_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC8_HPC,group.by = "group",label=T)




HCC9_HPC <- SCTransform(HCC9_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC9_HPC <- RunPCA(HCC9_HPC, verbose = FALSE)
HCC9_HPC <- RunUMAP(HCC9_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC9_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC9_HPC,group.by = "group",label=T)



sta_merge <- merge(HCC1_HPC, y = c(HCC3_HPC,HCC5_HPC,HCC6_HPC,HCC7_HPC))
sta_merge <-PrepSCTFindMarkers(sta_merge)
merge_sta_markers <- FindMarkers(sta_merge,ident.1 = "primary",ident.2 = "satellite",group.by = "satellite_region" )
merge_sta_markers$gene <- row.names(merge_sta_markers)
merge_sta_markers_sig <- subset(merge_sta_markers,subset = p_val_adj < 0.05)
merge_sta_markers <- mutate(merge_sta_markers,state = case_when(avg_log2FC>0&p_val_adj < 0.05 ~ "Up-regulated", # 上调
                                                                        avg_log2FC<0&p_val_adj < 0.05 ~ "Down-regulated", # 下调
                                                          TRUE ~ "Unchanged"))
merge_sta_markers <- mutate(merge_sta_markers,state2 = case_when(gene %in% pt_up_sum ~ "Up common", # 上调
                                                                TRUE ~ "Unchanged"))


merge_sta_markers_sig_up <- subset(merge_sta_markers_sig,subset = avg_log2FC > 0)
merge_sta_markers_sig_down <- subset(merge_sta_markers_sig,subset = avg_log2FC < 0)


merge_sta_up_go <- enrichGO(gene  = row.names(merge_sta_markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
merge_sta_up_go <- as.data.frame(merge_sta_up_go@result)
merge_sta_up_go [,"logp"] <- -log10(merge_sta_up_go$pvalue)
ggplot(data = merge_sta_up_go[1:20,])+
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
       title="merge_sta_up")

merge_sta_down_go <- enrichGO(gene  = row.names(merge_sta_markers_sig_down),
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
merge_sta_down_go <- as.data.frame(merge_sta_down_go@result)
merge_sta_down_go [,"logp"] <- -log10(merge_sta_down_go$pvalue)
ggplot(data = merge_sta_down_go[1:20,])+
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
       title="merge_sta_down")





sta_down_Venn <- list(hcc1 = row.names(hcc1_sta_markers_sig_down), hcc3 = row.names(hcc3_sta_markers_sig_down), hcc5 = row.names(hcc5_sta_markers_sig_down),
                    hcc6 = row.names(hcc6_sta_markers_sig_down),hcc7 = row.names(hcc7_sta_markers_sig_down))

sta_up_Venn <- list(hcc1 = row.names(hcc1_sta_markers_sig_up), hcc3 = row.names(hcc3_sta_markers_sig_up), hcc5 = row.names(hcc5_sta_markers_sig_up),
                      hcc6 = row.names(hcc6_sta_markers_sig_up),hcc7 = row.names(hcc7_sta_markers_sig_up))

venn.diagram(sta_down_Venn, filename = 'sta_down.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, 
             fill = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), cat.cex = 0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'),cex = 2.5 )

venn.diagram(sta_up_Venn, filename = 'sta_up.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, 
             fill = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), cat.cex =0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'),cex = 2.5 )


VlnPlot(HCC1_HPC,features = "VEGFA",group.by = "satellite_region")

pt_up_sum <- intersect(row.names(hcc1_sta_markers_sig_up),row.names(hcc6_sta_markers_sig_up))
pt_up_sum_go <- enrichGO(gene  =pt_up_sum,
                                OrgDb      = org.Hs.eg.db,
                                keyType    = 'SYMBOL',
                                ont        = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
pt_up_sum_go <- as.data.frame(pt_up_sum_go@result)
pt_up_sum_go [,"logp"] <- -log10(pt_up_sum_go$pvalue)
pt_up_sum_go[,11:12] <- as.numeric(str_split_fixed(pt_up_sum_go$GeneRatio,"/",2))
pt_up_sum_go$GeneRatio <- pt_up_sum_go[,11]/pt_up_sum_go[,12]
ggplot(data = pt_up_sum_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(GeneRatio),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="pt_up")


ggplot(merge_sta_markers, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(size = 0.4, aes(color = state)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue","grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="cr-GADD45A vs Vector diff gene")




hcc1_av <- AverageExpression(HCC1_HPC,group.by = "orig.ident",assays = 'RNA')
hcc1_av <- hcc1_av[[1]]
hcc1_cg=names(tail(sort(apply(hcc1_av, 1, sd)),1000))
pheatmap::pheatmap(cor(hcc1_av[hcc1_cg,],method = 'spearman'))




hcc2_av <- AverageExpression(HCC2_HPC,group.by = "orig.ident",assays = 'RNA')
hcc2_av <- hcc2_av[[1]]
hcc2_cg=names(tail(sort(apply(hcc2_av, 1, sd)),1000))
pheatmap::pheatmap(cor(hcc2_av[hcc2_cg,],method = 'spearman'))




hcc3_av <- AverageExpression(HCC3_HPC,group.by = "orig.ident",assays = 'RNA')
hcc3_av <- hcc3_av[[1]]
hcc3_cg=names(tail(sort(apply(hcc3_av, 1, sd)),1000))
pheatmap::pheatmap(cor(hcc3_av[hcc3_cg,],method = 'spearman'))





hcc1_all <- as.data.frame(GetAssayData(object = HCC1_HPC, slot = "data"))
hcc1_cor <- as.data.frame(cor(hcc1_all[,],method = 'pearson'))
hcc1_anno <- as.data.frame(HCC1_HPC$orig.ident)
hcc1_anno2 <- as.data.frame(HCC1_HPC$satellite_region)
colnames(hcc1_anno) <- "origin"
colnames(hcc1_anno2) <- "origin"
pheatmap(hcc1_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc1_anno,annotation_col = hcc1_anno,
         cluster_rows = F,cluster_cols = F)
pheatmap(hcc1_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc1_anno2,annotation_col = hcc1_anno2,
         cluster_rows = F,cluster_cols = F)
pheatmap(hcc1_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc1_anno2,annotation_col = hcc1_anno2,
         cluster_rows = T,cluster_cols = T)

table(HCC1_HPC$orig.ident)

hcc1_pt2_pt4_mean <- mean(unlist(hcc1_cor[149:1393,1394:2180]))
hcc1_pt2_pt5_mean <- mean(unlist(hcc1_cor[149:1393,2181:2522]))
hcc1_pt4_pt5_mean <- mean(unlist(hcc1_cor[1394:2180,2181:2522]))

hcc1_pt2_pt2_mean <- mean(unlist(hcc1_cor[149:1393,149:1393]))
hcc1_pt4_pt4_mean <- mean(unlist(hcc1_cor[1394:2180,1394:2180]))
hcc1_pt5_pt5_mean <- mean(unlist(hcc1_cor[2181:2522,2181:2522]))

hcc1_pt2_pt4 <- as.data.frame(unlist(hcc1_cor[149:1393,1394:2180]))
colnames(hcc1_pt2_pt4) <- "cor"
hcc1_pt2_pt4$pair <- "hcc1_pt2_pt4"
hcc1_pt2_pt4$type <- "cross"
hcc1_pt2_pt4$patient <- "hcc1"
hcc1_pt2_pt5 <- as.data.frame(unlist(hcc1_cor[149:1393,2181:2522]))
colnames(hcc1_pt2_pt5) <- "cor"
hcc1_pt2_pt5$pair <- "hcc1_pt2_pt5"
hcc1_pt2_pt5$type <- "cross"
hcc1_pt2_pt5$patient <- "hcc1"
hcc1_pt4_pt5 <- as.data.frame(unlist(hcc1_cor[1394:2180,2181:2522]))
colnames(hcc1_pt4_pt5) <- "cor"
hcc1_pt4_pt5$pair <- "hcc1_pt4_pt5"
hcc1_pt4_pt5$type <- "cross"
hcc1_pt4_pt5$patient <- "hcc1"

hcc1_pt2_pt2 <- as.data.frame(unlist(hcc1_cor[149:1393,149:1393]))
colnames(hcc1_pt2_pt2) <- "cor"
hcc1_pt2_pt2$pair <- "hcc1_pt2_pt2"
hcc1_pt2_pt2$type <- "self"
hcc1_pt2_pt2$patient <- "hcc1"
hcc1_pt4_pt4 <- as.data.frame(unlist(hcc1_cor[1394:2180,1394:2180]))
colnames(hcc1_pt4_pt4) <- "cor"
hcc1_pt4_pt4$pair <- "hcc1_pt4_pt4"
hcc1_pt4_pt4$type <- "self"
hcc1_pt4_pt4$patient <- "hcc1"
hcc1_pt5_pt5 <- as.data.frame(unlist(hcc1_cor[2181:2522,2181:2522]))
colnames(hcc1_pt5_pt5) <- "cor"
hcc1_pt5_pt5$pair <- "hcc1_pt5_pt5"
hcc1_pt5_pt5$type <- "self"
hcc1_pt5_pt5$patient <- "hcc1"

hcc1_cor_all <- as.data.frame(rbind(hcc1_pt2_pt4,hcc1_pt2_pt5,hcc1_pt4_pt5,
                                     hcc1_pt2_pt2,hcc1_pt4_pt4,hcc1_pt5_pt5))

ggplot(hcc1_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "t.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())




hcc2_all <- as.data.frame(GetAssayData(object = HCC2_HPC, slot = "data"))
hcc2_cor <- as.data.frame(cor(hcc2_all,method = 'pearson'))
hcc2_anno <- as.data.frame(HCC2_HPC$orig.ident)
colnames(hcc2_anno) <- "origin"
pheatmap(hcc2_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc2_anno,annotation_col = hcc2_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC2_HPC$orig.ident)

hcc2_pt1_pt2_mean <- mean(unlist(hcc2_cor[74:1723,1624:2871]))
hcc2_pt1_pt3_mean <- mean(unlist(hcc2_cor[74:1623,2872:3169]))
hcc2_pt1_pt4_mean <- mean(unlist(hcc2_cor[74:1623,3170:3255]))
hcc2_pt2_pt3_mean <- mean(unlist(hcc2_cor[1624:2871,2872:3169]))
hcc2_pt2_pt4_mean <- mean(unlist(hcc2_cor[1624:2871,3170:3255]))
hcc2_pt3_pt4_mean <- mean(unlist(hcc2_cor[2872:3169,3170:3255]))


hcc2_pt1_pt1_mean <- mean(unlist(hcc2_cor[74:1623,74:1623]))
hcc2_pt2_pt2_mean <- mean(unlist(hcc2_cor[1624:2871,1624:2871]))
hcc2_pt3_pt3_mean <- mean(unlist(hcc2_cor[2872:3169,2872:3169]))
hcc2_pt4_pt4_mean <- mean(unlist(hcc2_cor[3170:3255,3170:3255]))

hcc1_pt2_pt4_mean <- mean(unlist(hcc1_cor[149:1393,1394:2180]))
hcc1_pt2_pt5_mean <- mean(unlist(hcc1_cor[149:1393,2181:2522]))
hcc1_pt4_pt5_mean <- mean(unlist(hcc1_cor[1394:2180,2181:2522]))

hcc1_pt2_pt2_mean <- mean(unlist(hcc1_cor[149:1393,149:1393]))
hcc1_pt4_pt4_mean <- mean(unlist(hcc1_cor[1394:2180,1394:2180]))
hcc1_pt5_pt5_mean <- mean(unlist(hcc1_cor[2181:2522,2181:2522]))

hcc2_pt1_pt2_mean <- mean(unlist(hcc2_cor[74:1623,1624:2871]))
hcc2_pt1_pt3_mean <- mean(unlist(hcc2_cor[74:1623,2872:3169]))
hcc2_pt1_pt4_mean <- mean(unlist(hcc2_cor[74:1623,3170:3255]))
hcc2_pt2_pt3_mean <- mean(unlist(hcc2_cor[1624:2871,2872:3169]))
hcc2_pt2_pt4_mean <- mean(unlist(hcc2_cor[1624:2871,3170:3255]))
hcc2_pt3_pt4_mean <- mean(unlist(hcc2_cor[2872:3169,3170:3255]))


hcc2_pt1_pt1_mean <- mean(unlist(hcc2_cor[74:1623,74:1623]))
hcc2_pt2_pt2_mean <- mean(unlist(hcc2_cor[1624:2871,1624:2871]))
hcc2_pt3_pt3_mean <- mean(unlist(hcc2_cor[2872:3169,2872:3169]))
hcc2_pt4_pt4_mean <- mean(unlist(hcc2_cor[3170:3255,3170:3255]))



hcc2_pt1_pt2 <- as.data.frame(unlist(hcc2_cor[74:1723,1624:2871]))
colnames(hcc2_pt1_pt2) <- "cor"
hcc2_pt1_pt2$pair <- "hcc2_pt1_pt2"
hcc2_pt1_pt2$type <- "cross"
hcc2_pt1_pt2$patient <- "hcc2"

hcc2_pt1_pt3 <- as.data.frame(unlist(hcc2_cor[74:1623,2872:3169]))
colnames(hcc2_pt1_pt3) <- "cor"
hcc2_pt1_pt3$pair <- "hcc2_pt1_pt3"
hcc2_pt1_pt3$type <- "cross"
hcc2_pt1_pt3$patient <- "hcc2"

hcc2_pt1_pt4 <- as.data.frame(unlist(hcc2_cor[74:1623,3170:3255]))
colnames(hcc2_pt1_pt4) <- "cor"
hcc2_pt1_pt4$pair <- "hcc2_pt1_pt4"
hcc2_pt1_pt4$type <- "cross"
hcc2_pt1_pt4$patient <- "hcc2"

hcc2_pt2_pt3 <- as.data.frame(unlist(hcc2_cor[1624:2871,2872:3169]))
colnames(hcc2_pt2_pt3) <- "cor"
hcc2_pt2_pt3$pair <- "hcc2_pt2_pt3"
hcc2_pt2_pt3$type <- "cross"
hcc2_pt2_pt3$patient <- "hcc2"

hcc2_pt2_pt4 <- as.data.frame(unlist(hcc2_cor[1624:2871,3170:3255]))
colnames(hcc2_pt2_pt4) <- "cor"
hcc2_pt2_pt4$pair <- "hcc2_pt2_pt4"
hcc2_pt2_pt4$type <- "cross"
hcc2_pt2_pt4$patient <- "hcc2"

hcc2_pt3_pt4 <- as.data.frame(unlist(hcc2_cor[2872:3169,3170:3255]))
colnames(hcc2_pt3_pt4) <- "cor"
hcc2_pt3_pt4$pair <- "hcc2_pt3_pt4"
hcc2_pt3_pt4$type <- "cross"
hcc2_pt3_pt4$patient <- "hcc2"


hcc2_pt1_pt1 <- as.data.frame(unlist(hcc2_cor[74:1623,74:1623]))
colnames(hcc2_pt1_pt1) <- "cor"
hcc2_pt1_pt1$pair <- "hcc2_pt1_pt1"
hcc2_pt1_pt1$type <- "self"
hcc2_pt1_pt1$patient <- "hcc2"

hcc2_pt2_pt2 <- as.data.frame(unlist(hcc2_cor[1624:2871,1624:2871]))
colnames(hcc2_pt2_pt2) <- "cor"
hcc2_pt2_pt2$pair <- "hcc2_pt2_pt2"
hcc2_pt2_pt2$type <- "self"
hcc2_pt2_pt2$patient <- "hcc2"

hcc2_pt3_pt3 <- as.data.frame(unlist(hcc2_cor[2872:3169,2872:3169]))
colnames(hcc2_pt3_pt3) <- "cor"
hcc2_pt3_pt3$pair <- "hcc2_pt3_pt3"
hcc2_pt3_pt3$type <- "self"
hcc2_pt3_pt3$patient <- "hcc2"

hcc2_pt4_pt4 <- as.data.frame(unlist(hcc2_cor[3170:3255,3170:3255]))
colnames(hcc2_pt4_pt4) <- "cor"
hcc2_pt4_pt4$pair <- "hcc2_pt4_pt4"
hcc2_pt4_pt4$type <- "self"
hcc2_pt4_pt4$patient <- "hcc2"

hcc2_cor_all <- as.data.frame(rbind(hcc2_pt1_pt2,hcc2_pt1_pt3,hcc2_pt1_pt4,
                                    hcc2_pt2_pt3,hcc2_pt2_pt4,hcc2_pt3_pt4,
                                    hcc2_pt1_pt1,hcc2_pt2_pt2,hcc2_pt3_pt3,
                                    hcc2_pt4_pt4))

ggplot(hcc2_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "t.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc3_all <- as.data.frame(GetAssayData(object = HCC3_HPC, slot = "data"))
hcc3_cor <- as.data.frame(cor(hcc3_all[,],method = 'pearson'))
hcc3_anno <- as.data.frame(HCC3_HPC$orig.ident)
colnames(hcc3_anno) <- "origin"
pheatmap(hcc3_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc3_anno,annotation_col = hcc3_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC3_HPC$orig.ident)

hcc3_pt1_pt2_mean <- mean(unlist(hcc3_cor[13:156,157:245]))
hcc3_pt1_pt3_mean <- mean(unlist(hcc3_cor[13:156,245:248]))
hcc3_pt2_pt3_mean <- mean(unlist(hcc3_cor[157:245,245:248]))

hcc3_pt1_pt1_mean <- mean(unlist(hcc3_cor[13:156,13:156]))
hcc3_pt2_pt2_mean <- mean(unlist(hcc3_cor[157:245,157:245]))
hcc3_pt3_pt3_mean <- mean(unlist(hcc3_cor[245:248,245:248]))

hcc3_pt1_pt2 <- as.data.frame(unlist(hcc3_cor[13:156,157:245]))
colnames(hcc3_pt1_pt2) <- "cor"
hcc3_pt1_pt2$pair <- "hcc3_pt1_pt2"
hcc3_pt1_pt2$type <- "cross"
hcc3_pt1_pt2$patient <- "hcc3"
hcc3_pt1_pt3 <- as.data.frame(unlist(hcc3_cor[13:156,245:248]))
colnames(hcc3_pt1_pt3) <- "cor"
hcc3_pt1_pt3$pair <- "hcc3_pt1_pt3"
hcc3_pt1_pt3$type <- "cross"
hcc3_pt1_pt3$patient <- "hcc3"
hcc3_pt2_pt3 <- as.data.frame(unlist(hcc3_cor[157:245,245:248]))
colnames(hcc3_pt2_pt3) <- "cor"
hcc3_pt2_pt3$pair <- "hcc3_pt2_pt3"
hcc3_pt2_pt3$type <- "cross"
hcc3_pt2_pt3$patient <- "hcc3"

hcc3_pt1_pt1 <- as.data.frame(unlist(hcc3_cor[13:156,13:156]))
colnames(hcc3_pt1_pt1) <- "cor"
hcc3_pt1_pt1$pair <- "hcc3_pt1_pt1"
hcc3_pt1_pt1$type <- "self"
hcc3_pt1_pt1$patient <- "hcc3"
hcc3_pt2_pt2 <- as.data.frame(unlist(hcc3_cor[157:245,157:245]))
colnames(hcc3_pt2_pt2) <- "cor"
hcc3_pt2_pt2$pair <- "hcc3_pt2_pt2"
hcc3_pt2_pt2$type <- "self"
hcc3_pt2_pt2$patient <- "hcc3"
hcc3_pt3_pt3 <- as.data.frame(unlist(hcc3_cor[245:248,245:248]))
colnames(hcc3_pt3_pt3) <- "cor"
hcc3_pt3_pt3$pair <- "hcc3_pt3_pt3"
hcc3_pt3_pt3$type <- "self"
hcc3_pt3_pt3$patient <- "hcc3"

hcc3_cor_all <- as.data.frame(rbind(hcc3_pt1_pt2,hcc3_pt1_pt3,hcc3_pt2_pt3,
                                    hcc3_pt1_pt1,hcc3_pt2_pt2,hcc3_pt3_pt3))


ggplot(hcc3_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(hcc3_cor_all,aes(type,cor,color=type))+
  geom_violin(width=0.5)+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "t.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


hcc4_all <- as.data.frame(GetAssayData(object = HCC4_HPC, slot = "data"))
hcc4_cor <- as.data.frame(cor(hcc4_all[,],method = 'pearson'))
hcc4_anno <- as.data.frame(HCC4_HPC$orig.ident)
colnames(hcc4_anno) <- "origin"
pheatmap(hcc4_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc4_anno,annotation_col = hcc4_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC4_HPC$orig.ident)

hcc4_pt1_pt2_mean <- mean(unlist(hcc4_cor[3:16,17:61]))
hcc4_pt1_pt3_mean <- mean(unlist(hcc4_cor[3:16,62:69]))
hcc4_pt2_pt3_mean <- mean(unlist(hcc4_cor[17:61,62:69]))


hcc4_pt1_pt1_mean <- mean(unlist(hcc4_cor[3:16,3:16]))
hcc4_pt2_pt2_mean <- mean(unlist(hcc4_cor[17:61,17:61]))
hcc4_pt3_pt3_mean <- mean(unlist(hcc4_cor[62:69,62:69]))


hcc4_pt1_pt2 <- as.data.frame(unlist(hcc4_cor[3:16,17:61]))
colnames(hcc4_pt1_pt2) <- "cor"
hcc4_pt1_pt2$pair <- "hcc4_pt1_pt2"
hcc4_pt1_pt2$type <- "cross"
hcc4_pt1_pt2$patient <- "hcc4"
hcc4_pt1_pt3 <- as.data.frame(unlist(hcc4_cor[3:16,62:69]))
colnames(hcc4_pt1_pt3) <- "cor"
hcc4_pt1_pt3$pair <- "hcc4_pt1_pt3"
hcc4_pt1_pt3$type <- "cross"
hcc4_pt1_pt3$patient <- "hcc4"
hcc4_pt2_pt3 <- as.data.frame(unlist(hcc4_cor[17:61,62:69]))
colnames(hcc4_pt2_pt3) <- "cor"
hcc4_pt2_pt3$pair <- "hcc4_pt2_pt3"
hcc4_pt2_pt3$type <- "cross"
hcc4_pt2_pt3$patient <- "hcc4"

hcc4_pt1_pt1 <- as.data.frame(unlist(hcc4_cor[3:16,3:16]))
colnames(hcc4_pt1_pt1) <- "cor"
hcc4_pt1_pt1$pair <- "hcc4_pt1_pt1"
hcc4_pt1_pt1$type <- "self"
hcc4_pt1_pt1$patient <- "hcc4"
hcc4_pt2_pt2 <- as.data.frame(unlist(hcc4_cor[17:61,17:61]))
colnames(hcc4_pt2_pt2) <- "cor"
hcc4_pt2_pt2$pair <- "hcc4_pt2_pt2"
hcc4_pt2_pt2$type <- "self"
hcc4_pt2_pt2$patient <- "hcc4"
hcc4_pt3_pt3 <- as.data.frame(unlist(hcc4_cor[62:69,62:69]))
colnames(hcc4_pt3_pt3) <- "cor"
hcc4_pt3_pt3$pair <- "hcc4_pt3_pt3"
hcc4_pt3_pt3$type <- "self"
hcc4_pt3_pt3$patient <- "hcc4"


hcc4_cor_all <- as.data.frame(rbind(hcc4_pt1_pt2,hcc4_pt1_pt3,hcc4_pt2_pt3,
                                    hcc4_pt1_pt1,hcc4_pt2_pt2,hcc4_pt3_pt3))


ggplot(hcc4_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc5_all <- as.data.frame(GetAssayData(object = HCC5_HPC, slot = "data"))
hcc5_cor <- as.data.frame(cor(hcc5_all[,],method = 'pearson'))
hcc5_anno <- as.data.frame(HCC5_HPC$orig.ident)
colnames(hcc5_anno) <- "origin"
pheatmap(hcc5_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc5_anno,annotation_col = hcc5_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC5_HPC$orig.ident)

hcc5_pt1_pt2_mean <- mean(unlist(hcc5_cor[78:84,85:103]))
hcc5_pt1_pt3_mean <- mean(unlist(hcc5_cor[78:84,104:144]))
hcc5_pt1_pt4_mean <- mean(unlist(hcc5_cor[78:84,145:169]))
hcc5_pt1_pt5_mean <- mean(unlist(hcc5_cor[78:84,170:181]))
hcc5_pt2_pt3_mean <- mean(unlist(hcc5_cor[85:103,104:144]))
hcc5_pt2_pt4_mean <- mean(unlist(hcc5_cor[85:103,145:169]))
hcc5_pt2_pt5_mean <- mean(unlist(hcc5_cor[85:103,170:181]))
hcc5_pt3_pt4_mean <- mean(unlist(hcc5_cor[104:144,145:169]))
hcc5_pt3_pt5_mean <- mean(unlist(hcc5_cor[104:144,170:181]))
hcc5_pt4_pt5_mean <- mean(unlist(hcc5_cor[145:169,170:181]))


hcc5_pt1_pt1_mean <- mean(unlist(hcc5_cor[78:84,78:84]))
hcc5_pt2_pt2_mean <- mean(unlist(hcc5_cor[85:103,85:103]))
hcc5_pt3_pt3_mean <- mean(unlist(hcc5_cor[104:144,104:144]))
hcc5_pt4_pt4_mean <- mean(unlist(hcc5_cor[145:169,145:169]))
hcc5_pt5_pt5_mean <- mean(unlist(hcc5_cor[170:181,170:181]))



hcc5_pt1_pt2 <- as.data.frame(unlist(hcc5_cor[78:84,85:103]))
colnames(hcc5_pt1_pt2) <- "cor"
hcc5_pt1_pt2$pair <- "hcc5_pt1_pt2"
hcc5_pt1_pt2$type <- "cross"
hcc5_pt1_pt2$patient <- "hcc4"
hcc5_pt1_pt3 <- as.data.frame(unlist(hcc5_cor[78:84,104:144]))
colnames(hcc5_pt1_pt3) <- "cor"
hcc5_pt1_pt3$pair <- "hcc5_pt1_pt3"
hcc5_pt1_pt3$type <- "cross"
hcc5_pt1_pt3$patient <- "hcc5"
hcc5_pt1_pt4 <- as.data.frame(unlist(hcc5_cor[78:84,145:169]))
colnames(hcc5_pt1_pt4) <- "cor"
hcc5_pt1_pt4$pair <- "hcc5_pt1_pt4"
hcc5_pt1_pt4$type <- "cross"
hcc5_pt1_pt4$patient <- "hcc5"
hcc5_pt1_pt5 <- as.data.frame(unlist(hcc5_cor[78:84,170:181]))
colnames(hcc5_pt1_pt5) <- "cor"
hcc5_pt1_pt5$pair <- "hcc4_pt1_pt2"
hcc5_pt1_pt5$type <- "cross"
hcc5_pt1_pt5$patient <- "hcc5"
hcc5_pt2_pt3 <- as.data.frame(unlist(hcc5_cor[85:103,104:144]))
colnames(hcc5_pt2_pt3) <- "cor"
hcc5_pt2_pt3$pair <- "hcc5_pt2_pt3"
hcc5_pt2_pt3$type <- "cross"
hcc5_pt2_pt3$patient <- "hcc5"
hcc5_pt2_pt4 <- as.data.frame(unlist(hcc5_cor[85:103,145:169]))
colnames(hcc5_pt2_pt4) <- "cor"
hcc5_pt2_pt4$pair <- "hcc5_pt2_pt4"
hcc5_pt2_pt4$type <- "cross"
hcc5_pt2_pt4$patient <- "hcc5"
hcc5_pt2_pt5 <- as.data.frame(unlist(hcc5_cor[85:103,170:181]))
colnames(hcc5_pt2_pt5) <- "cor"
hcc5_pt2_pt5$pair <- "hcc5_pt2_pt5"
hcc5_pt2_pt5$type <- "cross"
hcc5_pt2_pt5$patient <- "hcc5"
hcc5_pt3_pt4 <- as.data.frame(unlist(hcc5_cor[104:144,145:169]))
colnames(hcc5_pt3_pt4) <- "cor"
hcc5_pt3_pt4$pair <- "hcc5_pt3_pt4" 
hcc5_pt3_pt4$type <- "cross"
hcc5_pt3_pt4$patient <- "hcc5"
hcc5_pt3_pt5 <- as.data.frame(unlist(hcc5_cor[104:144,170:181]))
colnames(hcc5_pt3_pt5) <- "cor"
hcc5_pt3_pt5$pair <- "hcc5_pt3_pt5"
hcc5_pt3_pt5$type <- "cross"
hcc5_pt3_pt5$patient <- "hcc5"
hcc5_pt4_pt5 <- as.data.frame(unlist(hcc5_cor[145:169,170:181]))
colnames(hcc5_pt4_pt5) <- "cor"
hcc5_pt4_pt5$pair <- "hcc5_pt4_pt5"
hcc5_pt4_pt5$type <- "cross"
hcc5_pt4_pt5$patient <- "hcc5"


hcc5_pt1_pt1 <- as.data.frame(unlist(hcc5_cor[78:84,78:84]))
colnames(hcc5_pt1_pt1) <- "cor"
hcc5_pt1_pt1$pair <- "hcc5_pt1_pt1"
hcc5_pt1_pt1$type <- "self"
hcc5_pt1_pt1$patient <- "hcc5"
hcc5_pt2_pt2 <- as.data.frame(unlist(hcc5_cor[85:103,85:103]))
colnames(hcc5_pt2_pt2) <- "cor"
hcc5_pt2_pt2$pair <- "hcc5_pt2_pt2"
hcc5_pt2_pt2$type <- "self"
hcc5_pt2_pt2$patient <- "hcc5"
hcc5_pt3_pt3 <- as.data.frame(unlist(hcc5_cor[104:144,104:144]))
colnames(hcc5_pt3_pt3) <- "cor"
hcc5_pt3_pt3$pair <- "hcc5_pt3_pt3"
hcc5_pt3_pt3$type <- "self"
hcc5_pt3_pt3$patient <- "hcc5"
hcc5_pt4_pt4 <- as.data.frame(unlist(hcc5_cor[145:169,145:169]))
colnames(hcc5_pt4_pt4) <- "cor"
hcc5_pt4_pt4$pair <- "hcc5_pt4_pt4"
hcc5_pt4_pt4$type <- "self"
hcc5_pt4_pt4$patient <- "hcc5"
hcc5_pt5_pt5 <- as.data.frame(unlist(hcc5_cor[170:181,170:181]))
colnames(hcc5_pt5_pt5) <- "cor"
hcc5_pt5_pt5$pair <- "hcc5_pt5_pt5"
hcc5_pt5_pt5$type <- "self"
hcc5_pt5_pt5$patient <- "hcc5"


hcc5_cor_all <- as.data.frame(rbind(hcc5_pt1_pt2, hcc5_pt1_pt3, hcc5_pt1_pt4, hcc5_pt1_pt5, hcc5_pt2_pt3, hcc5_pt2_pt4,
                                      hcc5_pt2_pt5, hcc5_pt3_pt4, hcc5_pt3_pt5, hcc5_pt4_pt5, hcc5_pt1_pt1, hcc5_pt2_pt2,
                                      hcc5_pt3_pt3, hcc5_pt4_pt4, hcc5_pt5_pt5))


ggplot(hcc5_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())



hcc6_all <- as.data.frame(GetAssayData(object = HCC6_HPC, slot = "data"))
hcc6_cor <- as.data.frame(cor(hcc6_all[,],method = 'pearson'))
hcc6_anno <- as.data.frame(HCC6_HPC$orig.ident)
colnames(hcc6_anno) <- "origin"
hcc6_anno2 <- as.data.frame(HCC6_HPC$satellite_region)
colnames(hcc6_anno2) <- "origin"
pheatmap(hcc6_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc6_anno,annotation_col = hcc6_anno,
         cluster_rows = F,cluster_cols = F)
table(HCC6_HPC$orig.ident)

hcc6_pt1_pt2_mean <- mean(unlist(hcc6_cor[19:115,116:169]))
hcc6_pt1_pt4_mean <- mean(unlist(hcc6_cor[19:115,170:298]))
hcc6_pt1_pt5_mean <- mean(unlist(hcc6_cor[19:115,299:429]))
hcc6_pt1_pt6_mean <- mean(unlist(hcc6_cor[19:115,430:434]))
hcc6_pt2_pt4_mean <- mean(unlist(hcc6_cor[116:169,170:298]))
hcc6_pt2_pt5_mean <- mean(unlist(hcc6_cor[116:169,299:429]))
hcc6_pt2_pt6_mean <- mean(unlist(hcc6_cor[116:169,104:144]))
hcc6_pt4_pt5_mean <- mean(unlist(hcc6_cor[104:144,170:298]))
hcc6_pt4_pt6_mean <- mean(unlist(hcc6_cor[104:144,430:434]))
hcc6_pt5_pt6_mean <- mean(unlist(hcc6_cor[170:298,430:434]))



hcc6_pt1_pt1_mean <- mean(unlist(hcc6_cor[19:115,19:115]))
hcc6_pt2_pt2_mean <- mean(unlist(hcc6_cor[116:169,116:169]))
hcc6_pt4_pt4_mean <- mean(unlist(hcc6_cor[104:144,104:144]))
hcc6_pt5_pt5_mean <- mean(unlist(hcc6_cor[170:298,170:298]))
hcc6_pt6_pt6_mean <- mean(unlist(hcc6_cor[430:434,430:434]))




hcc6_pt1_pt2 <- as.data.frame(unlist(hcc6_cor[19:115,116:169]))
colnames(hcc6_pt1_pt2) <- "cor"
hcc6_pt1_pt2$pair <- "hcc6_pt1_pt2"
hcc6_pt1_pt2$type <- "cross"
hcc6_pt1_pt2$patient <- "hcc6"
hcc6_pt1_pt4 <- as.data.frame(unlist(hcc6_cor[19:115,170:298]))
colnames(hcc6_pt1_pt4) <- "cor"
hcc6_pt1_pt4$pair <- "hcc6_pt1_pt4"
hcc6_pt1_pt4$type <- "cross"
hcc6_pt1_pt4$patient <- "hcc6"
hcc6_pt1_pt5 <- as.data.frame(unlist(hcc6_cor[19:115,299:429]))
colnames(hcc6_pt1_pt5) <- "cor"
hcc6_pt1_pt5$pair <- "hcc6_pt1_pt5"
hcc6_pt1_pt5$type <- "cross"
hcc6_pt1_pt5$patient <- "hcc6"
hcc6_pt1_pt6 <- as.data.frame(unlist(hcc6_cor[19:115,430:434]))
colnames(hcc6_pt1_pt6) <- "cor"
hcc6_pt1_pt6$pair <- "hcc6_pt1_pt6"
hcc6_pt1_pt6$type <- "cross"
hcc6_pt1_pt6$patient <- "hcc6"
hcc6_pt2_pt4 <- as.data.frame(unlist(hcc6_cor[116:169,170:298]))
colnames(hcc6_pt2_pt4) <- "cor"
hcc6_pt2_pt4$pair <- "hcc6_pt2_pt4"
hcc6_pt2_pt4$type <- "cross"
hcc6_pt2_pt4$patient <- "hcc6"
hcc6_pt2_pt5 <- as.data.frame(unlist(hcc6_cor[116:169,299:429]))
colnames(hcc6_pt2_pt5) <- "cor"
hcc6_pt2_pt5$pair <- "hcc6_pt2_pt5"
hcc6_pt2_pt5$type <- "cross"
hcc6_pt2_pt5$patient <- "hcc6"
hcc6_pt2_pt6 <- as.data.frame(unlist(hcc6_cor[116:169,104:144]))
colnames(hcc6_pt2_pt6) <- "cor"
hcc6_pt2_pt6$pair <- "hcc6_pt2_pt6"
hcc6_pt2_pt6$type <- "cross"
hcc6_pt2_pt6$patient <- "hcc6"
hcc6_pt4_pt5 <- as.data.frame(unlist(hcc6_cor[104:144,170:298]))
colnames(hcc6_pt4_pt5) <- "cor"
hcc6_pt4_pt5$pair <- "hcc6_pt4_pt5"
hcc6_pt4_pt5$type <- "cross"
hcc6_pt4_pt5$patient <- "hcc6"
hcc6_pt4_pt6 <- as.data.frame(unlist(hcc6_cor[104:144,430:434]))
colnames(hcc6_pt4_pt6) <- "cor"
hcc6_pt4_pt6$pair <- "hcc6_pt4_pt6"
hcc6_pt4_pt6$type <- "cross"
hcc6_pt4_pt6$patient <- "hcc6"
hcc6_pt5_pt6 <- as.data.frame(unlist(hcc6_cor[170:298,430:434]))
colnames(hcc6_pt5_pt6) <- "cor"
hcc6_pt5_pt6$pair <- "hcc6_pt5_pt6"
hcc6_pt5_pt6$type <- "cross"
hcc6_pt5_pt6$patient <- "hcc6"



hcc6_pt1_pt1 <- as.data.frame(unlist(hcc6_cor[19:115,19:115]))
colnames(hcc6_pt1_pt1) <- "cor"
hcc6_pt1_pt1$pair <- "hcc6_pt1_pt1"
hcc6_pt1_pt1$type <- "self"
hcc6_pt1_pt1$patient <- "hcc6"
hcc6_pt2_pt2 <- as.data.frame(unlist(hcc6_cor[116:169,116:169]))
colnames(hcc6_pt2_pt2) <- "cor"
hcc6_pt2_pt2$pair <- "hcc6_pt2_pt2"
hcc6_pt2_pt2$type <- "self"
hcc6_pt2_pt2$patient <- "hcc6"
hcc6_pt4_pt4 <- as.data.frame(unlist(hcc6_cor[104:144,104:144]))
colnames(hcc6_pt4_pt4) <- "cor"
hcc6_pt4_pt4$pair <- "hcc6_pt4_pt4"
hcc6_pt4_pt4$type <- "self"
hcc6_pt4_pt4$patient <- "hcc6"
hcc6_pt5_pt5 <- as.data.frame(unlist(hcc6_cor[170:298,170:298]))
colnames(hcc6_pt5_pt5) <- "cor"
hcc6_pt5_pt5$pair <- "hcc6_pt5_pt5"
hcc6_pt5_pt5$type <- "cross"
hcc6_pt5_pt5$patient <- "hcc6"
hcc6_pt6_pt6 <- as.data.frame(unlist(hcc6_cor[430:434,430:434]))
colnames(hcc6_pt6_pt6) <- "cor"
hcc6_pt6_pt6$pair <- "hcc6_pt6_pt6"
hcc6_pt6_pt6$type <- "self"
hcc6_pt6_pt6$patient <- "hcc6"

hcc6_cor_all <- as.data.frame(rbind(hcc6_pt1_pt2, hcc6_pt1_pt4, hcc6_pt1_pt5, hcc6_pt1_pt6, hcc6_pt2_pt4, hcc6_pt2_pt5,
                                    hcc6_pt2_pt6, hcc6_pt4_pt5, hcc6_pt4_pt6, hcc6_pt5_pt6, hcc6_pt1_pt1, hcc6_pt2_pt2,
                                    hcc6_pt4_pt4, hcc6_pt5_pt5, hcc6_pt6_pt6))





ggplot(hcc6_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


hcc7_all <- as.data.frame(GetAssayData(object = HCC7_HPC, slot = "data"))
hcc7_cor <- as.data.frame(cor(hcc7_all[,],method = 'pearson'))
hcc7_anno <- as.data.frame(HCC7_HPC$orig.ident)
colnames(hcc7_anno) <- "origin"
pheatmap(hcc7_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc7_anno,annotation_col = hcc7_anno,
         cluster_rows = F,cluster_cols = F)
table(HCC7_HPC$orig.ident)

hcc7_pt1_pt2_mean <- mean(unlist(hcc7_cor[11:42,43:58]))
hcc7_pt1_pt3_mean <- mean(unlist(hcc7_cor[11:42,59:116]))
hcc7_pt1_pt4_mean <- mean(unlist(hcc7_cor[11:42,117:159]))
hcc7_pt1_pt5_mean <- mean(unlist(hcc7_cor[11:42,160:174]))
hcc7_pt1_pt6_mean <- mean(unlist(hcc7_cor[11:42,175:187]))
hcc7_pt2_pt3_mean <- mean(unlist(hcc7_cor[43:58,59:116]))
hcc7_pt2_pt4_mean <- mean(unlist(hcc7_cor[43:58,117:159]))
hcc7_pt2_pt5_mean <- mean(unlist(hcc7_cor[43:58,160:174]))
hcc7_pt2_pt6_mean <- mean(unlist(hcc7_cor[43:58,175:187]))
hcc7_pt3_pt4_mean <- mean(unlist(hcc7_cor[59:116,117:159]))
hcc7_pt3_pt5_mean <- mean(unlist(hcc7_cor[59:116,160:174]))
hcc7_pt3_pt6_mean <- mean(unlist(hcc7_cor[59:116,175:187]))
hcc7_pt4_pt5_mean <- mean(unlist(hcc7_cor[117:159,160:174]))
hcc7_pt4_pt6_mean <- mean(unlist(hcc7_cor[117:159,175:187]))
hcc7_pt5_pt6_mean <- mean(unlist(hcc7_cor[160:174,175:187]))


hcc7_pt1_pt1_mean <- mean(unlist(hcc7_cor[11:42,11:42]))
hcc7_pt2_pt2_mean <- mean(unlist(hcc7_cor[43:58,43:58]))
hcc7_pt3_pt3_mean <- mean(unlist(hcc7_cor[59:116,59:116]))
hcc7_pt4_pt4_mean <- mean(unlist(hcc7_cor[117:159,117:159]))
hcc7_pt5_pt5_mean <- mean(unlist(hcc7_cor[160:174,160:174]))
hcc7_pt6_pt6_mean <- mean(unlist(hcc7_cor[175:187,175:187]))



hcc7_pt1_pt2 <- as.data.frame(unlist(hcc7_cor[11:42,43:58]))
colnames(hcc7_pt1_pt2) <- "cor"
hcc7_pt1_pt2$pair <- "hcc7_pt1_pt2"
hcc7_pt1_pt2$type <- "cross"
hcc7_pt1_pt2$patient <- "hcc7"
hcc7_pt1_pt3 <- as.data.frame(unlist(hcc7_cor[11:42,59:116]))
colnames(hcc7_pt1_pt3) <- "cor"
hcc7_pt1_pt3$pair <- "hcc7_pt1_pt3"
hcc7_pt1_pt3$type <- "cross"
hcc7_pt1_pt3$patient <- "hcc7"
hcc7_pt1_pt4 <- as.data.frame(unlist(hcc7_cor[11:42,117:159]))
colnames(hcc7_pt1_pt4) <- "cor"
hcc7_pt1_pt4$pair <- "hcc7_pt1_pt4"
hcc7_pt1_pt4$type <- "cross"
hcc7_pt1_pt4$patient <- "hcc7"
hcc7_pt1_pt5 <- as.data.frame(unlist(hcc7_cor[11:42,160:174]))
colnames(hcc7_pt1_pt5) <- "cor"
hcc7_pt1_pt5$pair <- "hcc7_pt1_pt5"
hcc7_pt1_pt5$type <- "cross"
hcc7_pt1_pt5$patient <- "hcc7"
hcc7_pt1_pt6 <- as.data.frame(unlist(hcc7_cor[11:42,175:187]))
colnames(hcc7_pt1_pt6) <- "cor"
hcc7_pt1_pt6$pair <- "hcc7_pt1_pt6"
hcc7_pt1_pt6$type <- "cross"
hcc7_pt1_pt6$patient <- "hcc7"
hcc7_pt2_pt3 <- as.data.frame(unlist(hcc7_cor[43:58,59:116]))
colnames(hcc7_pt2_pt3) <- "cor"
hcc7_pt2_pt3$pair <- "hcc7_pt2_pt3"
hcc7_pt2_pt3$type <- "cross"
hcc7_pt2_pt3$patient <- "hcc7"
hcc7_pt2_pt4 <- as.data.frame(unlist(hcc7_cor[43:58,117:159]))
colnames(hcc7_pt2_pt4) <- "cor"
hcc7_pt2_pt4$pair <- "hcc7_pt2_pt4"
hcc7_pt2_pt4$type <- "cross"
hcc7_pt2_pt4$patient <- "hcc7"
hcc7_pt2_pt5 <- as.data.frame(unlist(hcc7_cor[43:58,160:174]))
colnames(hcc7_pt2_pt5) <- "cor"
hcc7_pt2_pt5$pair <- "hcc7_pt2_pt5"
hcc7_pt2_pt5$type <- "cross"
hcc7_pt2_pt5$patient <- "hcc7"
hcc7_pt2_pt6 <- as.data.frame(unlist(hcc7_cor[43:58,175:187]))
colnames(hcc7_pt2_pt6) <- "cor"
hcc7_pt2_pt6$pair <- "hcc7_pt2_pt6"
hcc7_pt2_pt6$type <- "cross"
hcc7_pt2_pt6$patient <- "hcc7"
hcc7_pt3_pt4 <- as.data.frame(unlist(hcc7_cor[59:116,117:159]))
colnames(hcc7_pt3_pt4) <- "cor"
hcc7_pt3_pt4$pair <- "hcc7_pt3_pt4"
hcc7_pt3_pt4$type <- "cross"
hcc7_pt3_pt4$patient <- "hcc7"
hcc7_pt3_pt5 <- as.data.frame(unlist(hcc7_cor[59:116,160:174]))
colnames(hcc7_pt3_pt5) <- "cor"
hcc7_pt3_pt5$pair <- "hcc7_pt3_pt5"
hcc7_pt3_pt5$type <- "cross"
hcc7_pt3_pt5$patient <- "hcc7"
hcc7_pt3_pt6 <- as.data.frame(unlist(hcc7_cor[59:116,175:187]))
colnames(hcc7_pt3_pt6) <- "cor"
hcc7_pt3_pt6$pair <- "hcc7_pt3_pt6"
hcc7_pt3_pt6$type <- "cross"
hcc7_pt3_pt6$patient <- "hcc7"
hcc7_pt4_pt5 <- as.data.frame(unlist(hcc7_cor[117:159,160:174]))
colnames(hcc7_pt4_pt5) <- "cor"
hcc7_pt4_pt5$pair <- "hcc7_pt4_pt5"
hcc7_pt4_pt5$type <- "cross"
hcc7_pt4_pt5$patient <- "hcc7"
hcc7_pt4_pt6 <- as.data.frame(unlist(hcc7_cor[117:159,175:187]))
colnames(hcc7_pt4_pt6) <- "cor"
hcc7_pt4_pt6$pair <- "hcc7_pt4_pt6"
hcc7_pt4_pt6$type <- "cross"
hcc7_pt4_pt6$patient <- "hcc7"
hcc7_pt5_pt6 <- as.data.frame(unlist(hcc7_cor[160:174,175:187]))
colnames(hcc7_pt5_pt6) <- "cor"
hcc7_pt5_pt6$pair <- "hcc7_pt5_pt6"
hcc7_pt5_pt6$type <- "cross"
hcc7_pt5_pt6$patient <- "hcc7"



hcc7_pt1_pt1 <- as.data.frame(unlist(hcc7_cor[11:42,11:42]))
colnames(hcc7_pt1_pt1) <- "cor"
hcc7_pt1_pt1$pair <- "hcc7_pt1_pt1"
hcc7_pt1_pt1$type <- "self"
hcc7_pt1_pt1$patient <- "hcc7"
hcc7_pt2_pt2 <- as.data.frame(unlist(hcc7_cor[43:58,43:58]))
colnames(hcc7_pt2_pt2) <- "cor"
hcc7_pt2_pt2$pair <- "hcc7_pt2_pt2"
hcc7_pt2_pt2$type <- "self"
hcc7_pt2_pt2$patient <- "hcc7"
hcc7_pt3_pt3 <- as.data.frame(unlist(hcc7_cor[59:116,59:116]))
colnames(hcc7_pt3_pt3) <- "cor"
hcc7_pt3_pt3$pair <- "hcc7_pt3_pt3"
hcc7_pt3_pt3$type <- "self"
hcc7_pt3_pt3$patient <- "hcc7"
hcc7_pt4_pt4 <- as.data.frame(unlist(hcc7_cor[117:159,117:159]))
colnames(hcc7_pt4_pt4) <- "cor"
hcc7_pt4_pt4$pair <- "hcc7_pt4_pt4"
hcc7_pt4_pt4$type <- "self"
hcc7_pt4_pt4$patient <- "hcc7"
hcc7_pt5_pt5 <- as.data.frame(unlist(hcc7_cor[160:174,160:174]))
colnames(hcc7_pt5_pt5) <- "cor"
hcc7_pt5_pt5$pair <- "hcc7_pt5_pt5"
hcc7_pt5_pt5$type <- "self"
hcc7_pt5_pt5$patient <- "hcc7"
hcc7_pt6_pt6 <- as.data.frame(unlist(hcc7_cor[175:187,175:187]))
colnames(hcc7_pt6_pt6) <- "cor"
hcc7_pt6_pt6$pair <- "hcc7_pt6_pt6"
hcc7_pt6_pt6$type <- "self"
hcc7_pt6_pt6$patient <- "hcc7"


hcc7_cor_all <- as.data.frame(rbind(hcc7_pt1_pt2, hcc7_pt1_pt3, hcc7_pt1_pt4, hcc7_pt1_pt5, hcc7_pt1_pt6, hcc7_pt2_pt3,
                                    hcc7_pt2_pt4, hcc7_pt2_pt5, hcc7_pt2_pt6, hcc7_pt3_pt4, hcc7_pt3_pt5, hcc7_pt3_pt6,
                                    hcc7_pt4_pt5, hcc7_pt4_pt6, hcc7_pt5_pt6,
                                    hcc7_pt1_pt1, hcc7_pt2_pt2, hcc7_pt3_pt3, hcc7_pt4_pt4, hcc7_pt5_pt5, hcc7_pt6_pt6))


ggplot(hcc7_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc8_all <- as.data.frame(GetAssayData(object = HCC8_HPC, slot = "data"))
hcc8_cor <- as.data.frame(cor(hcc8_all[,],method = 'pearson'))
hcc8_anno <- as.data.frame(HCC8_HPC$orig.ident)
colnames(hcc8_anno) <- "origin"
pheatmap(hcc8_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc8_anno,annotation_col = hcc8_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC8_HPC$orig.ident)

hcc8_pt1_pt2_mean <- mean(unlist(hcc8_cor[78:1145,1146:2083]))
hcc8_pt1_pt4_mean <- mean(unlist(hcc8_cor[78:1145,2084:2114]))
hcc8_pt2_pt4_mean <- mean(unlist(hcc8_cor[1146:2083,2084:2114]))



hcc8_pt1_pt1_mean <- mean(unlist(hcc8_cor[78:1145,78:1145]))
hcc8_pt2_pt2_mean <- mean(unlist(hcc8_cor[1146:2083,1146:2083]))
hcc8_pt4_pt4_mean <- mean(unlist(hcc8_cor[2084:2114,2084:2114]))




hcc8_pt1_pt2 <- as.data.frame(unlist(hcc8_cor[78:1145,1146:2083]))
colnames(hcc8_pt1_pt2) <- "cor"
hcc8_pt1_pt2$pair <- "hcc8_pt1_pt2"
hcc8_pt1_pt2$type <- "cross"
hcc8_pt1_pt2$patient <- "hcc8"
hcc8_pt1_pt4 <- as.data.frame(unlist(hcc8_cor[78:1145,2084:2114]))
colnames(hcc8_pt1_pt4) <- "cor"
hcc8_pt1_pt4$pair <- "hcc8_pt1_pt4"
hcc8_pt1_pt4$type <- "cross"
hcc8_pt1_pt4$patient <- "hcc8"
hcc8_pt2_pt4 <- as.data.frame(unlist(hcc8_cor[1146:2083,2084:2114]))
colnames(hcc8_pt2_pt4) <- "cor"
hcc8_pt2_pt4$pair <- "hcc8_pt2_pt4"
hcc8_pt2_pt4$type <- "cross"
hcc8_pt2_pt4$patient <- "hcc8"


hcc8_pt1_pt1 <- as.data.frame(unlist(hcc8_cor[78:1145,78:1145]))
colnames(hcc8_pt1_pt1) <- "cor"
hcc8_pt1_pt1$pair <- "hcc8_pt1_pt1"
hcc8_pt1_pt1$type <- "self"
hcc8_pt1_pt1$patient <- "hcc8"
hcc8_pt2_pt2 <- as.data.frame(unlist(hcc8_cor[1146:2083,1146:2083]))
colnames(hcc8_pt2_pt2) <- "cor"
hcc8_pt2_pt2$pair <- "hcc8_pt2_pt2"
hcc8_pt2_pt2$type <- "self"
hcc8_pt2_pt2$patient <- "hcc8"
hcc8_pt4_pt4 <- as.data.frame(unlist(hcc8_cor[2084:2114,2084:2114]))
colnames(hcc8_pt4_pt4) <- "cor"
hcc8_pt4_pt4$pair <- "hcc8_pt4_pt4"
hcc8_pt4_pt4$type <- "self"
hcc8_pt4_pt4$patient <- "hcc8"




hcc8_cor_all <- as.data.frame(rbind(hcc8_pt1_pt2,hcc8_pt1_pt4,hcc8_pt2_pt4,
                                    hcc8_pt1_pt1,hcc8_pt2_pt2,hcc8_pt4_pt4))


ggplot(hcc8_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


hcc9_all <- as.data.frame(GetAssayData(object = HCC9_HPC, slot = "data"))
hcc9_cor <- as.data.frame(cor(hcc9_all[,],method = 'pearson'))
hcc9_anno <- as.data.frame(HCC9_HPC$orig.ident)
colnames(hcc9_anno) <- "origin"
pheatmap(hcc9_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc9_anno,annotation_col = hcc9_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC9_HPC$orig.ident)

hcc9_pt1_pt3_mean <- mean(unlist(hcc9_cor[164:494,495:733]))
hcc9_pt1_pt4_mean <- mean(unlist(hcc9_cor[164:494,734:1148]))
hcc9_pt3_pt4_mean <- mean(unlist(hcc9_cor[495:733,734:1148]))


hcc9_pt1_pt1_mean <- mean(unlist(hcc9_cor[164:494,164:494]))
hcc9_pt3_pt3_mean <- mean(unlist(hcc9_cor[495:733,495:733]))
hcc9_pt4_pt4_mean <- mean(unlist(hcc9_cor[734:1148,734:1148]))


hcc9_pt1_pt3 <- as.data.frame(unlist(hcc9_cor[164:494,495:733]))
colnames(hcc9_pt1_pt3) <- "cor"
hcc9_pt1_pt3$pair <- "hcc9_pt1_pt3"
hcc9_pt1_pt3$type <- "cross"
hcc9_pt1_pt3$patient <- "hcc9"
hcc9_pt1_pt4 <- as.data.frame(unlist(hcc9_cor[164:494,734:1148]))
colnames(hcc9_pt1_pt4) <- "cor"
hcc9_pt1_pt4$pair <- "hcc9_pt1_pt4"
hcc9_pt1_pt4$type <- "cross"
hcc9_pt1_pt4$patient <- "hcc9"
hcc9_pt3_pt4 <- as.data.frame(unlist(hcc9_cor[495:733,734:1148]))
colnames(hcc9_pt3_pt4) <- "cor"
hcc9_pt3_pt4$pair <- "hcc9_pt3_pt4"
hcc9_pt3_pt4$type <- "cross"
hcc9_pt3_pt4$patient <- "hcc9"


hcc9_pt1_pt1 <- as.data.frame(unlist(hcc9_cor[164:494,164:494]))
colnames(hcc9_pt1_pt1) <- "cor"
hcc9_pt1_pt1$pair <- "hcc9_pt1_pt1"
hcc9_pt1_pt1$type <- "self"
hcc9_pt1_pt1$patient <- "hcc9"
hcc9_pt3_pt3 <- as.data.frame(unlist(hcc9_cor[495:733,495:733]))
colnames(hcc9_pt3_pt3) <- "cor"
hcc9_pt3_pt3$pair <- "hcc9_pt3_pt3"
hcc9_pt3_pt3$type <- "self"
hcc9_pt3_pt3$patient <- "hcc9"
hcc9_pt4_pt4 <- as.data.frame(unlist(hcc9_cor[734:1148,734:1148]))
colnames(hcc9_pt4_pt4) <- "cor"
hcc9_pt4_pt4$pair <- "hcc9_pt4_pt4"
hcc9_pt4_pt4$type <- "self"
hcc9_pt4_pt4$patient <- "hcc9"




hcc9_cor_all <- as.data.frame(rbind(hcc9_pt1_pt3,hcc9_pt1_pt4,hcc9_pt3_pt4,
                                    hcc9_pt1_pt1,hcc9_pt3_pt3,hcc9_pt4_pt4))


ggplot(hcc9_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())



hcc1_cor_mean <- as.data.frame(rbind(hcc1_pt2_pt4_mean,hcc1_pt2_pt5_mean ,hcc1_pt4_pt5_mean,
                                     hcc1_pt2_pt2_mean,hcc1_pt4_pt4_mean,hcc1_pt5_pt5_mean))
colnames(hcc1_cor_mean) <- "cor"
hcc1_cor_mean$comb <- row.names(hcc1_cor_mean)
hcc1_cor_mean$sample <- "hcc1"
hcc1_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))





hcc2_cor_mean <- as.data.frame(rbind(hcc2_pt1_pt2_mean,hcc2_pt1_pt3_mean,hcc2_pt1_pt4_mean,
                                     hcc2_pt2_pt3_mean,hcc2_pt2_pt4_mean,hcc2_pt3_pt4_mean,
                                     hcc2_pt1_pt1_mean,hcc2_pt2_pt2_mean,hcc2_pt3_pt3_mean,
                                     hcc2_pt4_pt4_mean))
colnames(hcc2_cor_mean) <- "cor"
hcc2_cor_mean$comb <- row.names(hcc2_cor_mean)
hcc2_cor_mean$sample <- "hcc2"
hcc2_cor_mean$comb_type <- c(rep("cross",6),rep("self",4))





hcc3_cor_mean <- as.data.frame(rbind(hcc3_pt1_pt2_mean,hcc3_pt1_pt3_mean,hcc3_pt2_pt3_mean,
                                     hcc3_pt1_pt1_mean,hcc3_pt2_pt2_mean,hcc3_pt3_pt3_mean))
colnames(hcc3_cor_mean) <- "cor"
hcc3_cor_mean$comb <- row.names(hcc3_cor_mean)
hcc3_cor_mean$sample <- "hcc3"
hcc3_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))




hcc4_cor_mean <- as.data.frame(rbind(hcc4_pt1_pt2_mean,hcc4_pt1_pt3_mean,hcc4_pt2_pt3_mean,
                                     hcc4_pt1_pt1_mean,hcc4_pt2_pt2_mean,hcc4_pt3_pt3_mean))
colnames(hcc4_cor_mean) <- "cor"
hcc4_cor_mean$comb <- row.names(hcc4_cor_mean)
hcc4_cor_mean$sample <- "hcc4"
hcc4_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))






hcc5_cor_mean <- as.data.frame(rbind(hcc5_pt1_pt2_mean,hcc5_pt1_pt3_mean,hcc5_pt1_pt4_mean,hcc5_pt1_pt5_mean,
                                     hcc5_pt2_pt3_mean,hcc5_pt2_pt4_mean,hcc5_pt2_pt5_mean,hcc5_pt3_pt4_mean,
                                     hcc5_pt3_pt5_mean,hcc5_pt4_pt5_mean,
                                     hcc5_pt1_pt1_mean,hcc5_pt2_pt2_mean,hcc5_pt3_pt3_mean,hcc5_pt4_pt4_mean,
                                     hcc5_pt5_pt5_mean))
colnames(hcc5_cor_mean) <- "cor"
hcc5_cor_mean$comb <- row.names(hcc5_cor_mean)
hcc5_cor_mean$sample <- "hcc5"
hcc5_cor_mean$comb_type <- c(rep("cross",10),rep("self",5))





hcc6_cor_mean <- as.data.frame(rbind(hcc6_pt1_pt2_mean,hcc6_pt1_pt6_mean,hcc6_pt1_pt4_mean,hcc6_pt1_pt5_mean,
                                     hcc6_pt2_pt6_mean,hcc6_pt2_pt4_mean,hcc6_pt2_pt5_mean,hcc6_pt4_pt6_mean,
                                     hcc6_pt5_pt6_mean,hcc6_pt4_pt5_mean,
                                     hcc6_pt1_pt1_mean,hcc6_pt2_pt2_mean,hcc6_pt6_pt6_mean,hcc6_pt4_pt4_mean,
                                     hcc6_pt5_pt5_mean))
colnames(hcc6_cor_mean) <- "cor"
hcc6_cor_mean$comb <- row.names(hcc6_cor_mean)
hcc6_cor_mean$sample <- "hcc6"
hcc6_cor_mean$comb_type <- c(rep("cross",10),rep("self",5))



hcc7_cor_mean <- as.data.frame(rbind(hcc7_pt1_pt2_mean,hcc7_pt1_pt3_mean,hcc7_pt1_pt4_mean,hcc7_pt1_pt5_mean,
                                     hcc7_pt1_pt6_mean,hcc7_pt2_pt3_mean,hcc7_pt2_pt4_mean,hcc7_pt2_pt5_mean,
                                     hcc7_pt2_pt6_mean,hcc7_pt3_pt4_mean,hcc7_pt3_pt5_mean,hcc7_pt3_pt6_mean,
                                     hcc7_pt4_pt5_mean,hcc7_pt4_pt6_mean,hcc7_pt5_pt6_mean,
                                     hcc7_pt1_pt1_mean,hcc7_pt2_pt2_mean,hcc7_pt3_pt3_mean,hcc7_pt4_pt4_mean,
                                     hcc7_pt5_pt5_mean,hcc7_pt6_pt6_mean))
colnames(hcc7_cor_mean) <- "cor"
hcc7_cor_mean$comb <- row.names(hcc7_cor_mean)
hcc7_cor_mean$sample <- "hcc7"
hcc7_cor_mean$comb_type <- c(rep("cross",15),rep("self",6))


hcc8_cor_mean <- as.data.frame(rbind(hcc8_pt1_pt2_mean,hcc8_pt1_pt4_mean,hcc8_pt2_pt4_mean,
                                     hcc8_pt1_pt1_mean,hcc8_pt2_pt2_mean,hcc8_pt4_pt4_mean))
colnames(hcc8_cor_mean) <- "cor"
hcc8_cor_mean$comb <- row.names(hcc8_cor_mean)
hcc8_cor_mean$sample <- "hcc8"
hcc8_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))


hcc9_cor_mean <- as.data.frame(rbind(hcc9_pt1_pt4_mean,hcc9_pt1_pt3_mean,hcc9_pt3_pt4_mean,
                                     hcc9_pt1_pt1_mean,hcc9_pt4_pt4_mean,hcc9_pt3_pt3_mean))
colnames(hcc9_cor_mean) <- "cor"
hcc9_cor_mean$comb <- row.names(hcc9_cor_mean)
hcc9_cor_mean$sample <- "hcc9"
hcc9_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))




sc_cor_sum <- rbind(hcc1_cor_mean,hcc2_cor_mean,hcc3_cor_mean,
                    hcc4_cor_mean,hcc5_cor_mean,hcc6_cor_mean,
                    hcc7_cor_mean,hcc8_cor_mean,hcc9_cor_mean)


ggplot(sc_cor_sum,aes(x=sample,y=cor,fill=comb_type,color=comb_type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")


VlnPlot(HPC,features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,split.by = "lib.method",group.by = "patient")






ggplot(sc_cor_sum,aes(comb_type,cor,color=comb_type))+
  geom_boxplot(width=0.5)+
  geom_point()+
  geom_jitter(width = 0.25)+
  facet_grid(~sample)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "t.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


sc_cor_all <- rbind(hcc1_cor_all,hcc2_cor_all,hcc3_cor_all,
                    hcc4_cor_all,hcc5_cor_all,hcc6_cor_all,
                    hcc7_cor_all,hcc8_cor_all,hcc9_cor_all)

ggplot(sc_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(sc_cor_all,aes(type,cor,color=type))+
  geom_violin(width=0.5,outlier.size=0)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc1_cross_mean <- mean(subset(hcc1_cor_all,subset=type=="cross")$cor)
hcc1_self_mean <- mean(subset(hcc1_cor_all,subset=type=="self")$cor)
hcc2_cross_mean <- mean(subset(hcc2_cor_all,subset=type=="cross")$cor)
hcc2_self_mean <- mean(subset(hcc2_cor_all,subset=type=="self")$cor)
hcc3_cross_mean <- mean(subset(hcc3_cor_all,subset=type=="cross")$cor)
hcc3_self_mean <- mean(subset(hcc3_cor_all,subset=type=="self")$cor)
hcc4_cross_mean <- mean(subset(hcc4_cor_all,subset=type=="cross")$cor)
hcc4_self_mean <- mean(subset(hcc4_cor_all,subset=type=="self")$cor)
hcc5_cross_mean <- mean(subset(hcc5_cor_all,subset=type=="cross")$cor)
hcc5_self_mean <- mean(subset(hcc5_cor_all,subset=type=="self")$cor)
hcc6_cross_mean <- mean(subset(hcc6_cor_all,subset=type=="cross")$cor)
hcc6_self_mean <- mean(subset(hcc6_cor_all,subset=type=="self")$cor)
hcc7_cross_mean <- mean(subset(hcc7_cor_all,subset=type=="cross")$cor)
hcc7_self_mean <- mean(subset(hcc7_cor_all,subset=type=="self")$cor)
hcc8_cross_mean <- mean(subset(hcc8_cor_all,subset=type=="cross")$cor)
hcc8_self_mean <- mean(subset(hcc8_cor_all,subset=type=="self")$cor)
hcc9_cross_mean <- mean(subset(hcc9_cor_all,subset=type=="cross")$cor)
hcc9_self_mean <- mean(subset(hcc9_cor_all,subset=type=="self")$cor)

hcc1_mean <- as.data.frame(abs(hcc1_self_mean-hcc1_cross_mean))
colnames(hcc1_mean) <- "cor_delta"
hcc1_mean$patient <- "hcc1"
hcc1_mean$type <- "CMN"


hcc2_mean <- as.data.frame(abs(hcc2_self_mean-hcc2_cross_mean))
colnames(hcc2_mean) <- "cor_delta"
hcc2_mean$patient <- "hcc2"
hcc2_mean$type <- "CMN"

hcc3_mean <- as.data.frame(abs(hcc3_self_mean-hcc3_cross_mean))
colnames(hcc3_mean) <- "cor_delta"
hcc3_mean$patient <- "hcc3"
hcc3_mean$type <- "SN"


hcc4_mean <- as.data.frame(abs(hcc4_self_mean-hcc4_cross_mean))
colnames(hcc4_mean) <- "cor_delta"
hcc4_mean$patient <- "hcc4"
hcc4_mean$type <- "CMN"


hcc5_mean <- as.data.frame(abs(hcc5_self_mean-hcc5_cross_mean))
colnames(hcc5_mean) <- "cor_delta"
hcc5_mean$patient <- "hcc5"
hcc5_mean$type <- "SN"

hcc6_mean <- as.data.frame(abs(hcc6_self_mean-hcc6_cross_mean))
colnames(hcc6_mean) <- "cor_delta"
hcc6_mean$patient <- "hcc6"
hcc6_mean$type <- "CMN"


hcc7_mean <- as.data.frame(abs(hcc7_self_mean-hcc7_cross_mean))
colnames(hcc7_mean) <- "cor_delta"
hcc7_mean$patient <- "hcc7"
hcc7_mean$type <- "SN"


hcc8_mean <- as.data.frame(abs(hcc8_self_mean-hcc8_cross_mean))
colnames(hcc8_mean) <- "cor_delta"
hcc8_mean$patient <- "hcc8"
hcc8_mean$type <- "CMN"



hcc9_mean <- as.data.frame(abs(hcc9_self_mean-hcc9_cross_mean))
colnames(hcc9_mean) <- "cor_delta"
hcc9_mean$patient <- "hcc9"
hcc9_mean$type <- "CMN"

sc_mean_sum <- rbind(hcc1_mean,hcc2_mean,hcc3_mean,
                    hcc4_mean,hcc5_mean,hcc6_mean,
                    hcc7_mean,hcc8_mean,hcc9_mean)




ggplot(sc_mean_sum,aes(type,cor_delta,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  geom_jitter(width=0.25)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("CMN","SN")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =0.2)+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


######satellite#####

######hcc1,hcc3,hcc5,hcc6,hcc7####
######hcc1-pt2/pt4####pt5#####

hcc1_primary <- subset(HCC1_HPC,subset=satellite_region=="primary")
hcc1_satellite <- subset(HCC1_HPC,subset=satellite_region=="satellite")
hcc1_primary_mtx <- as.data.frame(GetAssayData(object = hcc1_primary, slot = "data"))

hcc1_primary_cor <- as.data.frame(cor(hcc1_primary_mtx,method = 'pearson'))
hcc1_satellite_mtx <- as.data.frame(GetAssayData(object = hcc1_satellite, slot = "data"))
hcc1_satellite_cor <- as.data.frame(cor(hcc1_satellite_mtx,method = 'pearson'))

pheatmap(hcc1_satellite_cor)

hcc1_primary_cor2 <- as.data.frame(unlist(unique(hcc1_primary_cor)))
colnames(hcc1_primary_cor2) <- "cor"
hcc1_primary_cor2$type <- "primary"
hcc1_satellite_cor2 <- as.data.frame(unlist(unique(hcc1_satellite_cor)))
colnames(hcc1_satellite_cor2) <- "cor"
hcc1_satellite_cor2$type <- "satellite"


hcc1_pt2 <- subset(HCC1_HPC, subset = orig.ident=="PT2")
hcc1_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt2, slot = "data"))
hcc1_pt2_top_prm_genes <- names(tail(sort(apply(hcc1_pt2_mtx,1,sd)),3000))
hcc1_pt2_primary_top_cor <- as.data.frame(cor(hcc1_pt2_mtx[hcc1_pt2_top_prm_genes,],method = 'pearson'))
hcc1_pt2_cor <- as.data.frame(cor(hcc1_pt2_mtx,method = 'pearson'))

for(i in 1:nrow(hcc1_pt2_primary_top_cor)){
  vector <- as.vector(hcc1_pt2_primary_top_cor[i,1:(nrow(hcc1_pt2_primary_top_cor)-i+1)])
  if(i == 1){
    hcc1_pt2_cor_nu <- vector
  }
  else{
    hcc1_pt2_cor_nu <- c(hcc1_pt2_cor_nu,vector)
  }
}



hcc1_pt4 <- subset(HCC1_HPC, subset = orig.ident=="PT4")
hcc1_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt4, slot = "data"))
hcc1_pt4_top_prm_genes <- names(tail(sort(apply(hcc1_pt4_mtx,1,sd)),3000))
hcc1_pt4_primary_top_cor <- as.data.frame(cor(hcc1_pt4_mtx[hcc1_pt4_top_prm_genes,],method = 'pearson'))
hcc1_pt4_cor <- as.data.frame(cor(hcc1_pt4_mtx,method = 'pearson'))
for(i in 1:nrow(hcc1_pt4_primary_top_cor)){
  vector <- as.vector(hcc1_pt4_primary_top_cor[i,1:(nrow(hcc1_pt4_primary_top_cor)-i+1)])
  if(i == 1){
    hcc1_pt4_cor_nu <- vector
  }
  else{
    hcc1_pt4_cor_nu <- c(hcc1_pt4_cor_nu,vector)
  }
}

hcc1_pt5 <- subset(HCC1_HPC, subset = orig.ident=="PT5")
hcc1_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt5, slot = "data"))
hcc1_pt5_top_prm_genes <- names(tail(sort(apply(hcc1_pt5_mtx,1,sd)),3000))
hcc1_pt5_primary_top_cor <- as.data.frame(cor(hcc1_pt5_mtx[hcc1_pt4_top_prm_genes,],method = 'pearson'))
hcc1_pt5_cor <- as.data.frame(cor(hcc1_pt5_mtx,method = 'pearson'))
for(i in 1:nrow(hcc1_pt5_primary_top_cor)){
  vector <- as.vector(hcc1_pt5_primary_top_cor[i,1:(nrow(hcc1_pt5_primary_top_cor)-i+1)])
  if(i == 1){
    hcc1_pt5_cor_nu <- vector
  }
  else{
    hcc1_pt5_cor_nu <- c(hcc1_pt5_cor_nu,vector)
  }
}


hcc1_primary_cor3 <- as.data.frame(c(hcc1_pt4_cor_nu,hcc1_pt2_cor_nu))
hcc1_primary_cor3 <-as.data.frame(t(hcc1_primary_cor3))
colnames(hcc1_primary_cor3) <- "cor"
hcc1_primary_cor3$type <- "primary"
hcc1_satellite_cor3 <- as.data.frame(unlist(unique(hcc1_pt5_cor_nu)))
colnames(hcc1_satellite_cor3) <- "cor"
hcc1_satellite_cor3$type <- "satellite"
hcc1_ps_info2 <- rbind(hcc1_primary_cor3,hcc1_satellite_cor3)


hcc1_ps_info <- rbind(hcc1_primary_cor2,hcc1_satellite_cor2)
ggplot(data=hcc1_ps_info2,aes(x=cor,stat(density),color=type))+
  geom_freqpoly(binwidth=0.03,linewidth=1)+
  theme_classic()+
  scale_x_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,3))+
  labs(y="Frequency",x="Cor")+
  scale_color_manual(values = c("#2d2884","#c2a20c"),
                     name="Element")

ggplot(data=hcc1_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC1")





hcc1_avg_sta <- AverageExpression(HCC1_HPC,
                                  group.by = "satellite_region",
                                  assays = "RNA")
hcc1_avg_sta <- hcc1_avg_sta[[1]]
head(hcc1_avg_sta)
hcc1_top_sta_genes <- names(tail(sort(apply(hcc1_avg_sta,1,sd)),3000))
hcc1_top_sta_cor <- cor(hcc1_avg_sta[hcc1_top_sta_genes,],method = 'spearman')
pheatmap(hcc1_top_sta_cor)

hcc1_primary_mtx <- as.data.frame(GetAssayData(object = hcc1_primary, slot = "data"))
hcc1_top_prm_genes <- names(tail(sort(apply(hcc1_primary_mtx,1,sd)),3000))
hcc1_primary_top_cor <- as.data.frame(cor(hcc1_primary_mtx[hcc1_top_prm_genes,],method = 'pearson'))
hcc1_satellite_mtx <- as.data.frame(GetAssayData(object = hcc1_satellite, slot = "data"))
hcc1_top_sta_genes <- names(tail(sort(apply(hcc1_satellite_mtx,1,sd)),3000))
hcc1_satellite_top_cor <- as.data.frame(cor(hcc1_satellite_mtx[hcc1_top_prm_genes,],method = 'pearson'))

hcc1_primary_top_cor2 <- as.data.frame(unlist(unique(hcc1_primary_top_cor)))
colnames(hcc1_primary_top_cor2) <- "cor"
hcc1_primary_top_cor2$type <- "primary"
hcc1_satellite_top_cor2 <- as.data.frame(unlist(unique(hcc1_satellite_top_cor)))
colnames(hcc1_satellite_top_cor2) <- "cor"
hcc1_satellite_top_cor2$type <- "satellite"

hcc1_ps_top_info <- rbind(hcc1_primary_top_cor2,hcc1_satellite_top_cor2)
ggplot(data=hcc1_ps_top_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC1")


######hcc3-pt1/pt2####pt3#####
hcc3_primary <- subset(HCC3_HPC,subset=satellite_region=="primary")
hcc3_satellite <- subset(HCC3_HPC,subset=satellite_region=="satellite")
hcc3_primary_mtx <- as.data.frame(GetAssayData(object = hcc3_primary, slot = "data"))
hcc3_primary_cor <- as.data.frame(cor(hcc3_primary_mtx,method = 'pearson'))
hcc3_satellite_mtx <- as.data.frame(GetAssayData(object = hcc3_satellite, slot = "data"))
hcc3_satellite_cor <- as.data.frame(cor(hcc3_satellite_mtx,method = 'pearson'))

hcc3_primary_cor2 <- as.data.frame(unlist(unique(hcc3_primary_cor)))
colnames(hcc3_primary_cor2) <- "cor"
hcc3_primary_cor2$type <- "primary"
hcc3_satellite_cor2 <- as.data.frame(unlist(unique(hcc3_satellite_cor)))
colnames(hcc3_satellite_cor2) <- "cor"
hcc3_satellite_cor2$type <- "satellite"

hcc3_ps_info <- rbind(hcc3_primary_cor2,hcc3_satellite_cor2)
ggplot(data=hcc3_ps_info,aes(x=cor,stat(density),color=type))+
  geom_freqpoly(binwidth=0.03,linewidth=1)+
  theme_classic()+
  scale_x_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,3))+
  labs(y="Frequency",x="Cor",title = "HCC3")+
  scale_color_manual(values = c("#2d2884","#c2a20c"),
                     name="Element")

ggplot(data=hcc3_ps_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC3")



hcc3_primary_mtx <- as.data.frame(GetAssayData(object = hcc3_primary, slot = "data"))
hcc3_top_prm_genes <- names(tail(sort(apply(hcc3_primary_mtx,1,sd)),3000))
hcc3_primary_top_cor <- as.data.frame(cor(hcc3_primary_mtx[hcc3_top_prm_genes,],method = 'pearson'))
hcc3_satellite_mtx <- as.data.frame(GetAssayData(object = hcc3_satellite, slot = "data"))
hcc3_top_sta_genes <- names(tail(sort(apply(hcc3_satellite_mtx,1,sd)),3000))
hcc3_satellite_top_cor <- as.data.frame(cor(hcc3_satellite_mtx[hcc3_top_prm_genes,],method = 'pearson'))

hcc3_primary_top_cor2 <- as.data.frame(unlist(unique(hcc3_primary_top_cor)))
colnames(hcc3_primary_top_cor2) <- "cor"
hcc3_primary_top_cor2$type <- "primary"
hcc3_satellite_top_cor2 <- as.data.frame(unlist(unique(hcc3_satellite_top_cor)))
colnames(hcc3_satellite_top_cor2) <- "cor"
hcc3_satellite_top_cor2$type <- "satellite"

hcc3_ps_top_info <- rbind(hcc3_primary_top_cor2,hcc3_satellite_top_cor2)
ggplot(data=hcc3_ps_top_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC3")



hcc3_pt1 <- subset(HCC3_HPC, subset = orig.ident=="PT1")
hcc3_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt1, slot = "data"))
hcc3_pt1_cor <- as.data.frame(cor(hcc3_pt1_mtx,method = 'pearson'))
hcc3_pt1_cor_nu <- unique(as.vector(as.matrix(hcc3_pt1_cor)))
for(i in 1:nrow(hcc3_pt1_cor)){
  vector <- as.vector(hcc3_pt1_cor[i,1:(nrow(hcc3_pt1_cor)-i+1)])
  if(i == 1){
    hcc3_pt1_cor_nu <- vector
  }
  else{
    hcc3_pt1_cor_nu <- c(hcc3_pt1_cor_nu,vector)
  }
}

hcc3_pt2 <- subset(HCC3_HPC, subset = orig.ident=="PT2")
hcc3_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt2, slot = "data"))
hcc3_pt2_cor <- as.data.frame(cor(hcc3_pt2_mtx,method = 'pearson'))
hcc3_pt2_cor_nu <- unique(as.vector(as.matrix(hcc3_pt2_cor)))
for(i in 1:nrow(hcc3_pt2_cor)){
  vector <- as.vector(hcc3_pt2_cor[i,1:(nrow(hcc3_pt2_cor)-i+1)])
  if(i == 1){
    hcc3_pt2_cor_nu <- vector
  }
  else{
    hcc3_pt2_cor_nu <- c(hcc3_pt2_cor_nu,vector)
  }
}

hcc3_pt3 <- subset(HCC3_HPC, subset = orig.ident=="PT3")
hcc3_pt3_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt3, slot = "data"))
hcc3_pt3_cor <- as.data.frame(cor(hcc3_pt3_mtx,method = 'pearson'))
hcc3_pt3_cor_nu <- unique(as.vector(as.matrix(hcc3_pt3_cor)))
for(i in 1:nrow(hcc3_pt3_cor)){
  vector <- as.vector(hcc3_pt3_cor[i,1:(nrow(hcc3_pt3_cor)-i+1)])
  if(i == 1){
    hcc3_pt3_cor_nu <- vector
  }
  else{
    hcc3_pt3_cor_nu <- c(hcc3_pt3_cor_nu,vector)
  }
}

hcc3_primary_cor3 <- as.data.frame(c(hcc3_pt1_cor_nu,hcc3_pt2_cor_nu))
hcc3_primary_cor3 <-as.data.frame(t(hcc3_primary_cor3))
colnames(hcc3_primary_cor3) <- "cor"
hcc3_primary_cor3$type <- "primary"
hcc3_satellite_cor3 <- as.data.frame(unlist(unique(hcc3_pt3_cor_nu)))
colnames(hcc3_satellite_cor3) <- "cor"
hcc3_satellite_cor3$type <- "satellite"
hcc3_ps_info2 <- rbind(hcc3_primary_cor3,hcc3_satellite_cor3)
ggplot(data=hcc3_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC3")



######hcc5-pt1/pt2/pt3/pt4####pt5#####
hcc5_primary <- subset(HCC5_HPC,subset=satellite_region=="primary")
hcc5_satellite <- subset(HCC5_HPC,subset=satellite_region=="satellite")
hcc5_primary_mtx <- as.data.frame(GetAssayData(object = hcc5_primary, slot = "data"))
hcc5_primary_cor <- as.data.frame(cor(hcc5_primary_mtx,method = 'pearson'))
hcc5_satellite_mtx <- as.data.frame(GetAssayData(object = hcc5_satellite, slot = "data"))
hcc5_satellite_cor <- as.data.frame(cor(hcc5_satellite_mtx,method = 'pearson'))

hcc5_primary_cor2 <- as.data.frame(unlist(unique(hcc5_primary_cor)))
colnames(hcc5_primary_cor2) <- "cor"
hcc5_primary_cor2$type <- "primary"
hcc5_satellite_cor2 <- as.data.frame(unlist(unique(hcc5_satellite_cor)))
colnames(hcc5_satellite_cor2) <- "cor"
hcc5_satellite_cor2$type <- "satellite"

hcc5_ps_info <- rbind(hcc5_primary_cor2,hcc5_satellite_cor2)
ggplot(data=hcc5_ps_info,aes(x=cor,stat(density),color=type))+
  geom_freqpoly(binwidth=0.03,linewidth=1)+
  theme_classic()+
  scale_x_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,3))+
  labs(y="Frequency",x="Cor")+
  scale_color_manual(values = c("#2d2884","#c2a20c"),
                     name="Element")

ggplot(data=hcc5_ps_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC5")




hcc5_primary_mtx <- as.data.frame(GetAssayData(object = hcc5_primary, slot = "data"))
hcc5_top_prm_genes <- names(tail(sort(apply(hcc5_primary_mtx,1,sd)),3000))
hcc5_primary_top_cor <- as.data.frame(cor(hcc5_primary_mtx[hcc5_top_prm_genes,],method = 'pearson'))
hcc5_satellite_mtx <- as.data.frame(GetAssayData(object = hcc5_satellite, slot = "data"))
hcc5_top_sta_genes <- names(tail(sort(apply(hcc5_satellite_mtx,1,sd)),3000))
hcc5_satellite_top_cor <- as.data.frame(cor(hcc5_satellite_mtx[hcc5_top_prm_genes,],method = 'pearson'))

hcc5_primary_top_cor2 <- as.data.frame(unlist(unique(hcc5_primary_top_cor)))
colnames(hcc5_primary_top_cor2) <- "cor"
hcc5_primary_top_cor2$type <- "primary"
hcc5_satellite_top_cor2 <- as.data.frame(unlist(unique(hcc5_satellite_top_cor)))
colnames(hcc5_satellite_top_cor2) <- "cor"
hcc5_satellite_top_cor2$type <- "satellite"

hcc5_ps_top_info <- rbind(hcc5_primary_top_cor2,hcc5_satellite_top_cor2)
ggplot(data=hcc5_ps_top_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC5")
######hcc6-pt1/pt2####pt4/pt5/pt6#####
hcc6_primary <- subset(HCC6_HPC,subset=satellite_region=="primary")
hcc6_satellite <- subset(HCC6_HPC,subset=satellite_region=="satellite")
hcc6_primary_mtx <- as.data.frame(GetAssayData(object = hcc6_primary, slot = "data"))
hcc6_primary_cor <- as.data.frame(cor(hcc6_primary_mtx,method = 'pearson'))
hcc6_satellite_mtx <- as.data.frame(GetAssayData(object = hcc6_satellite, slot = "data"))
hcc6_satellite_cor <- as.data.frame(cor(hcc6_satellite_mtx,method = 'pearson'))

hcc6_primary_cor2 <- as.data.frame(unlist(unique(hcc6_primary_cor)))
colnames(hcc6_primary_cor2) <- "cor"
hcc6_primary_cor2$type <- "primary"
hcc6_satellite_cor2 <- as.data.frame(unlist(unique(hcc6_satellite_cor)))
colnames(hcc6_satellite_cor2) <- "cor"
hcc6_satellite_cor2$type <- "satellite"

hcc6_ps_info <- rbind(hcc6_primary_cor2,hcc6_satellite_cor2)
ggplot(data=hcc6_ps_info,aes(x=cor,stat(density),color=type))+
  geom_freqpoly(binwidth=0.03,linewidth=1)+
  theme_classic()+
  scale_x_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,3))+
  labs(y="Frequency",x="Cor")+
  scale_color_manual(values = c("#2d2884","#c2a20c"),
                     name="Element")

ggplot(data=hcc6_ps_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC6")





hcc6_primary_mtx <- as.data.frame(GetAssayData(object = hcc6_primary, slot = "data"))
hcc6_top_prm_genes <- names(tail(sort(apply(hcc6_primary_mtx,1,sd)),3000))
hcc6_primary_top_cor <- as.data.frame(cor(hcc6_primary_mtx[hcc6_top_prm_genes,],method = 'pearson'))
hcc6_satellite_mtx <- as.data.frame(GetAssayData(object = hcc6_satellite, slot = "data"))
hcc6_top_sta_genes <- names(tail(sort(apply(hcc6_satellite_mtx,1,sd)),3000))
hcc6_satellite_top_cor <- as.data.frame(cor(hcc6_satellite_mtx[hcc6_top_prm_genes,],method = 'pearson'))

hcc6_primary_top_cor2 <- as.data.frame(unlist(unique(hcc6_primary_top_cor)))
colnames(hcc6_primary_top_cor2) <- "cor"
hcc6_primary_top_cor2$type <- "primary"
hcc6_satellite_top_cor2 <- as.data.frame(unlist(unique(hcc6_satellite_top_cor)))
colnames(hcc6_satellite_top_cor2) <- "cor"
hcc6_satellite_top_cor2$type <- "satellite"

hcc6_ps_top_info <- rbind(hcc6_primary_top_cor2,hcc6_satellite_top_cor2)
ggplot(data=hcc6_ps_top_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC6")

######hcc7-pt1/pt2/pt3/pt4/pt5###pt6##
hcc7_primary <- subset(HCC7_HPC,subset=satellite_region=="primary")
hcc7_satellite <- subset(HCC7_HPC,subset=satellite_region=="satellite")
hcc7_primary_mtx <- as.data.frame(GetAssayData(object = hcc7_primary, slot = "data"))
hcc7_primary_cor <- as.data.frame(cor(hcc7_primary_mtx,method = 'pearson'))
hcc7_satellite_mtx <- as.data.frame(GetAssayData(object = hcc7_satellite, slot = "data"))
hcc7_satellite_cor <- as.data.frame(cor(hcc7_satellite_mtx,method = 'pearson'))

hcc7_primary_cor2 <- as.data.frame(unlist(unique(hcc7_primary_cor)))
colnames(hcc7_primary_cor2) <- "cor"
hcc7_primary_cor2$type <- "primary"
hcc7_satellite_cor2 <- as.data.frame(unlist(unique(hcc7_satellite_cor)))
colnames(hcc7_satellite_cor2) <- "cor"
hcc7_satellite_cor2$type <- "satellite"

hcc7_ps_info <- rbind(hcc7_primary_cor2,hcc7_satellite_cor2)
ggplot(data=hcc7_ps_info,aes(x=cor,stat(density),color=type))+
  geom_freqpoly(binwidth=0.03,linewidth=1)+
  theme_classic()+
  scale_x_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,3))+
  labs(y="Frequency",x="Cor")+
  scale_color_manual(values = c("#2d2884","#c2a20c"),
                     name="Element")

ggplot(data=hcc7_ps_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC7")


hcc7_primary_mtx <- as.data.frame(GetAssayData(object = hcc7_primary, slot = "data"))
hcc7_top_prm_genes <- names(tail(sort(apply(hcc7_primary_mtx,1,sd)),3000))
hcc7_primary_top_cor <- as.data.frame(cor(hcc7_primary_mtx[hcc7_top_prm_genes,],method = 'pearson'))
hcc7_satellite_mtx <- as.data.frame(GetAssayData(object = hcc7_satellite, slot = "data"))
hcc7_top_sta_genes <- names(tail(sort(apply(hcc7_satellite_mtx,1,sd)),3000))
hcc7_satellite_top_cor <- as.data.frame(cor(hcc7_satellite_mtx[hcc7_top_prm_genes,],method = 'pearson'))

hcc7_primary_top_cor2 <- as.data.frame(unlist(unique(hcc7_primary_top_cor)))
colnames(hcc7_primary_top_cor2) <- "cor"
hcc7_primary_top_cor2$type <- "primary"
hcc7_satellite_top_cor2 <- as.data.frame(unlist(unique(hcc7_satellite_top_cor)))
colnames(hcc7_satellite_top_cor2) <- "cor"
hcc7_satellite_top_cor2$type <- "satellite"

hcc7_ps_top_info <- rbind(hcc7_primary_top_cor2,hcc7_satellite_top_cor2)
ggplot(data=hcc7_ps_top_info,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC7")
 




hpc <- as.data.frame(rbind(c(2522,2032,342),c(248,232,4),c(181,92,12),c(434,151,265),c(187,163,13)))
colnames(hpc) <- c("normal","primary","satellite")
rownames(hpc) <- c("HCC1","HCC3","HCC5","HCC6","HCC7")





tx_HPC <- subset(HPC,subset = lib.method =="10x")
tx_HPC <- subset(tx_HPC, subset = patient != "HCC1")
tx_HPC <- SCTransform(tx_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
tx_HPC <- RunPCA(tx_HPC, verbose = FALSE)
tx_HPC <- RunUMAP(tx_HPC, dims = 1:30, verbose = FALSE)

DimPlot(tx_HPC,group.by = "patient")
DimPlot(tx_HPC,group.by = "sample")
DimPlot(tx_HPC,group.by = "group")

tx_HPC$sample <- paste(tx_HPC$patient,tx_HPC$orig.ident,sep = "_")
tx_HPC_sample <- tx_HPC$sample
tx_HPC_methinfo = case_when(
  tx_HPC_sample  %in% c("HCC2_PT1","HCC2_PT2","HCC8_PT1","HCC8_PT2","HCC8_PT4","HCC9_PT1","HCC9_PT3","HCC9_PT4")~"Hypo-Methy",
  tx_HPC_sample  %in% c("HCC2_PT3","HCC2_PT4","HCC2_NT","HCC8_NT","HCC9_NT")~"Hyper-Methy",
  TRUE ~ as.character(tx_HPC_sample))
tx_HPC$meth_type <- tx_HPC_methinfo

DimPlot(tx_HPC,group.by = "meth_type")



immune_merge <- subset(bigseu, subset = new.ident %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                              "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                              "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                              "PlasmaB cell","Neutrophil"))


immune_merge_cmn <- subset(immune_merge, subset = cntype %in% "CMN")
immune_merge_sn <- subset(immune_merge, subset = cntype %in% "SN")

cell.prop_sn <- as.data.frame(prop.table(table(immune_merge_sn$new.ident,immune_merge_sn$cntype)))
cell.prop_sn <-  subset(cell.prop_sn, Var1 %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                  "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                  "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                  "PlasmaB cell","Neutrophil"))
cell.prop_cmn <- as.data.frame(prop.table(table(immune_merge_cmn$new.ident,immune_merge_cmn$cntype)))
cell.prop_cmn <-  subset(cell.prop_cmn, Var1 %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                  "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                  "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                  "PlasmaB cell","Neutrophil"))

cell.prop_cob <- rbind(select(cell.prop_sn,1,3,2),select(cell.prop_cmn,1,3,2))
colnames(cell.prop_cob) <- c("cell_type","Freq","CN_type")

ggplot(cell.prop_cob)+
  geom_bar(aes(x=CN_type,y=Freq,fill=CN_type),stat="identity", position=position_dodge())+
  facet_wrap(~ cell_type)


immune_merge$patient_pt <- paste(immune_merge$patient,immune_merge$orig.ident,sep = "_")
HPC$patient_pt <- paste(HPC$patient,HPC$orig.ident,sep = "_")
bigseu$patient_pt <- paste(bigseu$patient,bigseu$orig.ident,sep = "_")

immune_merge_nt <- subset(immune_merge, subset = group %in% "NT")
immune_merge_pt <- subset(immune_merge, subset = group %in% "P")
immune_merge_st <- subset(immune_merge, subset = group %in% "S")


cell.prop_nt <- as.data.frame(prop.table(table(immune_merge_nt$new.ident,immune_merge_nt$group)))
cell.prop_nt2 <- as.data.frame(prop.table(table(immune_merge_nt$new.ident,immune_merge_nt$patient_pt)))
cell.prop_nt <-  subset(cell.prop_nt, Var1 %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                  "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                  "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                  "PlasmaB cell","Neutrophil"))
cell.prop_pt <- as.data.frame(prop.table(table(immune_merge_pt$new.ident,immune_merge_pt$group)))
cell.prop_pt <-  subset(cell.prop_pt, Var1 %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                    "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                    "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                    "PlasmaB cell","Neutrophil"))
cell.prop_st <- as.data.frame(prop.table(table(immune_merge_st$new.ident,immune_merge_st$group)))
cell.prop_st <-  subset(cell.prop_st, Var1 %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                  "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                  "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                  "PlasmaB cell","Neutrophil"))

cell.prop_stpt_cob <- rbind(select(cell.prop_st,1,3,2),select(cell.prop_pt,1,3,2))
colnames(cell.prop_stpt_cob) <- c("cell_type","Freq","type")
ggplot(cell.prop_stpt_cob)+
  geom_bar(aes(x=type,y=Freq,fill=type),stat="identity", position=position_dodge())+
  facet_wrap(~ cell_type)

cell.prop_ntpt_cob <- rbind(select(cell.prop_nt,1,3,2),select(cell.prop_pt,1,3,2))
colnames(cell.prop_ntpt_cob) <- c("cell_type","Freq","type")

ggplot(cell.prop_ntpt_cob)+
  geom_bar(aes(x=type,y=Freq,fill=type),stat="identity", position=position_dodge())+
  facet_wrap(~ cell_type)

T_list <-names(table(subset(bigseu,subset = group != "NT")$patient_pt))
imtm_table <- as.data.frame(cbind(table(immune_merge$patient_pt),table(HPC$patient_pt)))

imtm_table <- imtm_table[T_list,]

colnames(imtm_table) <- c("immune","tumor")
imtm_table$radio <- imtm_table[,1]/imtm_table[,2]
imtm_table$patient_pt <- rownames(imtm_table)
imtm_table$patient <- str_split_fixed(imtm_table$patient_pt,"_",2)[,1]

ggplot(imtm_table,aes(x=patient,y=radio,fill=patient))+
  geom_boxplot(outlier.size = 0)+
  ylab('immune/tumor_ratio')+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=1)+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12,angle=45,color="black"),axis.text.y = element_text(size = 10,color="black"))

imtm_table2 <- fread("immune_ratio.txt")
imtm_table2[,'ratio'] <- log(imtm_table2$`cd45+`/imtm_table2$`cd45-`)

ggplot(imtm_table2,aes(x=Patient,y=ratio,fill=Patient))+
  geom_boxplot(outlier.size = 0)+
  ylab('log(immune/tumor_ratio)')+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=1)+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12,angle=45,color="black"),axis.text.y = element_text(size = 10,color="black"))

immune_merge_nt_result <- data.frame()

# 查找所有样品来源
nt_samples <- unique(immune_merge_nt@meta.data$patient_pt)

# 对每个样品来源，计算每种细胞类型的比例
for(s in nt_samples){
  # 选取特定的样品来源
  cells_of_sample <- subset(immune_merge_nt@meta.data, patient_pt == s)
  
  # 计算百分比
  prop <- prop.table(table(cells_of_sample$new.ident))
  
  # 将结果转换为数据框，并添加样品来源列
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'nt',
    'patient' = unique(cells_of_sample$patient)
  )
  
  # 将结果添加到最终结果数据框
  immune_merge_nt_result <- rbind(immune_merge_nt_result, prop_df)
}


immune_merge_pt_result <- data.frame()
pt_samples <- unique(immune_merge_pt@meta.data$patient_pt)

for(s in pt_samples){
  # 选取特定的样品来源
  cells_of_sample <- subset(immune_merge_pt@meta.data, patient_pt == s)
  
  # 计算百分比
  prop <- prop.table(table(cells_of_sample$new.ident))
  
  # 将结果转换为数据框，并添加样品来源列
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'pt',
    'patient' = unique(cells_of_sample$patient)
  )
  
  # 将结果添加到最终结果数据框
  immune_merge_pt_result <- rbind(immune_merge_pt_result, prop_df)
}



immune_merge_st_result <- data.frame()
st_samples <- unique(immune_merge_st@meta.data$patient_pt)

for(s in st_samples){
  # 选取特定的样品来源
  cells_of_sample <- subset(immune_merge_st@meta.data, patient_pt == s)
  
  # 计算百分比
  prop <- prop.table(table(cells_of_sample$new.ident))
  
  # 将结果转换为数据框，并添加样品来源列
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'st',
    'patient' = unique(cells_of_sample$patient)
  )
  
  # 将结果添加到最终结果数据框
  immune_merge_st_result <- rbind(immune_merge_st_result, prop_df)
}


prop_result_ptst_merge <- rbind(immune_merge_pt_result,immune_merge_st_result)
prop_result_ptnt_merge <- rbind(immune_merge_pt_result,immune_merge_nt_result)

prop_result_ptst_merge <- subset(prop_result_ptst_merge, subset = celltype %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                         "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                         "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                         "PlasmaB cell","Neutrophil"))
prop_result_ptnt_merge <- subset(prop_result_ptnt_merge, subset = celltype %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                                                  "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                                                  "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                                                  "PlasmaB cell","Neutrophil"))
prop_result_ptst_merge$celltype <-factor(prop_result_ptst_merge$celltype, levels = c( "CD8+ exhausted", "CD4+ memory", "CD4+ Treg","CD8+ cytotoxic",
                                                                                   "B cell", "CD8+ memory", "Dendritic cell","NK",
                                                                               "Neutrophil","Macrophage", "PlasmaB cell","Proliferative T"))

ggplot(prop_result_ptst_merge,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.2)+
  stat_compare_means(comparisons = list(c("pt","st")),
                     method = "t.test",label = "p.signif",
                     label.y = 0.13 )





prop_result_ptnt_merge_rm89  <- subset(prop_result_ptnt_merge,subset= patient %in% c('HCC1','HCC2','HCC3','HCC4','HCC5','HCC6','HCC7'))
ggplot(prop_result_ptnt_merge,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.size = 0.1)+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.4)+
  stat_compare_means(comparisons = list(c("pt","nt")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 0.32)


ggplot(prop_result_ptnt_merge_rm89,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.size = 0.1)+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.2)+
  stat_compare_means(comparisons = list(c("pt","nt")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 0.15 )


immune_merge_cmn2 <- subset(immune_merge_cmn, subset = group != "NT")
immune_merge_sn2 <- subset(immune_merge_sn, subset = group != "NT")


immune_merge_cmn_result <- data.frame()
cmn_samples <- unique(immune_merge_cmn@meta.data$patient_pt)

for(s in cmn_samples){
  cells_of_sample <- subset(immune_merge_cmn@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'cmn'
  )

  immune_merge_cmn_result <- rbind(immune_merge_cmn_result, prop_df)
}


immune_merge_sn_result <- data.frame()
sn_samples <- unique(immune_merge_sn@meta.data$patient_pt)

for(s in sn_samples){
  cells_of_sample <- subset(immune_merge_sn@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'sn'
  )
  
  immune_merge_sn_result <- rbind(immune_merge_sn_result, prop_df)
}



immune_prop_cmnsn_result_merge <- rbind(immune_merge_sn_result,immune_merge_cmn_result)


immune_prop_cmnsn_result_merge <- subset(immune_prop_cmnsn_result_merge, subset = celltype %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                                               "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                                               "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                                               "PlasmaB cell","Neutrophil"))

ggplot(immune_prop_cmnsn_result_merge,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.5)+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y = 0.4 )


immune_merge_cmn2 <- subset(immune_merge_cmn, subset = group != "NT")
immune_merge_sn2 <- subset(immune_merge_sn, subset = group != "NT")

immune_merge_cmn_result2 <- data.frame()
cmn_samples2 <- unique(immune_merge_cmn2@meta.data$patient_pt)

for(s in cmn_samples2){
  cells_of_sample <- subset(immune_merge_cmn2@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'cmn'
  )
  
  immune_merge_cmn_result2 <- rbind(immune_merge_cmn_result2, prop_df)
}


immune_merge_sn_result2 <- data.frame()
sn_samples2 <- unique(immune_merge_sn2@meta.data$patient_pt)

for(s in sn_samples2){
  cells_of_sample <- subset(immune_merge_sn2@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'sn'
  )
  
  immune_merge_sn_result2 <- rbind(immune_merge_sn_result2, prop_df)
}



immune_prop_cmnsn_result_merge2 <- rbind(immune_merge_sn_result2,immune_merge_cmn_result2)


immune_prop_cmnsn_result_merge2 <- subset(immune_prop_cmnsn_result_merge2, subset = celltype %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                                                     "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                                                     "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                                                     "PlasmaB cell","Neutrophil"))



immune_prop_cmnsn_result_merge2$celltype <-factor(immune_prop_cmnsn_result_merge2$celltype, levels = c( "CD8+ exhausted", "CD4+ memory", "CD4+ Treg","CD8+ cytotoxic",
                                                                                      "B cell", "CD8+ memory", "Dendritic cell","NK",
                                                                                      "Neutrophil","Macrophage", "PlasmaB cell","Proliferative T"))





ggplot(immune_prop_cmnsn_result_merge2,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.2)+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y = 0.13)









hpc_markers <- FindMarkers(HPC,ident.1 = "NT",ident.2 = "P",group.by = "group" )
hpc_markers_sig <- subset(hpc_markers,subset = p_val_adj < 0.05)
hpc_markers_sig_up <- subset(hpc_markers_sig,subset = avg_log2FC > 0)
hpc_markers_sig_down <- subset(hpc_markers_sig,subset = avg_log2FC < 0)

hpc_markers_up_go <- enrichGO(gene  = row.names(hpc_markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
hpc_markers_up_go <- as.data.frame(hpc_markers_up_go@result)
hpc_markers_up_go [,"logp"] <- -log10(hpc_markers_up_go$pvalue)
hpc_markers_up_go[,11:12] <- as.numeric(str_split_fixed(hpc_markers_up_go$GeneRatio,"/",2))
hpc_markers_up_go$GeneRatio <- hpc_markers_up_go[,11]/hpc_markers_up_go[,12]
ggplot(data = hpc_markers_up_go[1:20,])+
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
       title="pt_up")
ggplot(data = hpc_markers_up_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(GeneRatio),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="pt_up")



hpc_markers_down_go <- enrichGO(gene  = row.names(hpc_markers_sig_down),
                              OrgDb      = org.Hs.eg.db,
                              keyType    = 'SYMBOL',
                              ont        = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
hpc_markers_down_go <- as.data.frame(hpc_markers_down_go@result)
hpc_markers_down_go [,"logp"] <- -log10(hpc_markers_down_go$pvalue)
hpc_markers_down_go[,11:12] <- as.numeric(str_split_fixed(hpc_markers_down_go$GeneRatio,"/",2))
hpc_markers_down_go$GeneRatio <- hpc_markers_down_go[,11]/hpc_markers_down_go[,12]
ggplot(data = hpc_markers_down_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(GeneRatio),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="pt_down")
















immune_merge_hcc89 <- subset(immune_merge, subset = patient %in% c("HCC8","HCC9"))
immune_merge_hcc89 <- subset(immune_merge_hcc89,subset = group != "NT")
cell.prop_hcc89 <- as.data.frame(prop.table(table(immune_merge_hcc89$patient_pt,immune_merge_hcc89$new.ident)))

colnames(cell.prop_hcc89) <- c("HCC","origin","proportion")
ggplot(cell.prop_hcc89,aes(HCC,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  scale_fill_manual(values=my36colors)+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 14,color="black",angle = 45,hjust = 1),axis.text.y = element_text(size = 10,color="black"))









cell.prop_immune_merge <- as.data.frame(prop.table(table(immune_merge$patient_pt,immune_merge$new.ident)))

colnames(cell.prop_immune_merge) <- c("HCC","origin","proportion")

cell.prop_immune_merge <- subset(cell.prop_immune_merge, subset = origin %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                                                "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                                                "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                                                "PlasmaB cell","Neutrophil"))



cell.prop_immune_merge[,'patient'] <- str_split_fixed(cell.prop_immune_merge$HCC,"_",2)[,1]
cell.prop_immune_merge[,'sample'] <- str_split_fixed(cell.prop_immune_merge$HCC,"_",2)[,2]

hcc1_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC1')
hcc2_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC2')
hcc3_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC3')
hcc4_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC4')
hcc5_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC5')
hcc6_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC6')
hcc7_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC7')
hcc8_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC8')
hcc9_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC9')

immune_merge_data_list <- list(hcc1_cell.prop_immune_merge, hcc2_cell.prop_immune_merge, hcc3_cell.prop_immune_merge,
                               hcc4_cell.prop_immune_merge, hcc5_cell.prop_immune_merge, hcc6_cell.prop_immune_merge,
                               hcc7_cell.prop_immune_merge, hcc8_cell.prop_immune_merge, hcc9_cell.prop_immune_merge)



immuneprop_plot <- function(data) {
  ggplot(data,aes(sample,proportion,fill=origin))+
    geom_bar(stat="identity",position="fill")+
    ggtitle("")+
    scale_fill_manual(values=my36colors)+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'))+
    guides(fill=guide_legend(title=NULL))+
    labs(x=paste(data[1,4]))+
    theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 0),legend.text=element_text(size = 12))+
    theme(axis.text.x = element_text(size = 10,color="black",angle = 45,hjust = 1),axis.text.y = element_text(size = 10,color="black"))+
    NoLegend()
}

immuneprop_plot_list <- lapply(immune_merge_data_list, immuneprop_plot)

do.call(grid.arrange, c(immuneprop_plot_list, ncol = 5))









ggplot(cell.prop_immune_merge,aes(sample,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  scale_fill_manual(values=my36colors)+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 0),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10,color="black",angle = 45,hjust = 1),axis.text.y = element_text(size = 10,color="black"))




ggplot(hcc9_cell.prop_immune_merge,aes(sample,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  scale_fill_manual(values=my36colors)+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  labs(x=paste(hcc9_cell.prop_immune_merge[1,4]))+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 0),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10,color="black",angle = 45,hjust = 1),axis.text.y = element_text(size = 10,color="black"))
