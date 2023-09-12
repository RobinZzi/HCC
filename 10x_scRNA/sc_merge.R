rm(list = ls())
setwd("~/projects/hcc/analysis/merged_scrna")

save.image("merged_scRNA.Rdata")


bigseu <- readRDS("~/projects/hcc/data/10x_scRNA/merge/HCC-scRNA-seq-inte.rds")

bigseu$lib.method
DimPlot(bigseu,group.by = "group") #NT,PT,ST
DimPlot(bigseu,group.by = "lib.method")#10x,dropseq

DimPlot(bigseu,group.by = "orig.ident")# PT1/PT2/PT3

DimPlot(bigseu,group.by = "new.ident")#cell_type
DimPlot(bigseu,group.by = "patient")#patient

HPC <- subset(bigseu,subset=new.ident=="HPC")

HPC$patient

DimPlot(HPC,group.by = "patient")

hpc_merge <- SCTransform(HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
hpc_merge <- RunPCA(hpc_merge, verbose = FALSE)
hpc_merge <- RunUMAP(hpc_merge, dims = 1:30, verbose = FALSE)

DimPlot(hpc_merge,group.by = "patient",label=T)
DimPlot(hpc_merge,group.by = "group",label=T)
DimPlot(hpc_merge,group.by = "lib.method",label=T)



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



markers <- c('CD3D','CCL5', 'NKG7', 'GZMA', 'IL32', 'CD4','CD8A','FOXP3', 'CD3E', 'LTB', 'S100A8', 'S100A9', 'CD79A', "ENG" ,
             'FCN1', 'MS4A1', 'SPON2','FCER1A','SERPINF1', "LYZ","AMBP","COL1A2")
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置

VlnPlot(bigseu, features = markers,pt.size=0, stack = T,
        cols = my36colors,group.by = "new.ident")+NoLegend()+
  theme(axis.text.x = element_blank())

element_text(angle = 45,vjust = 0.5,hjust = 0.5)

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




HCC2_HPC <- SCTransform(HCC2_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC2_HPC <- RunPCA(HCC2_HPC, verbose = FALSE)
HCC2_HPC <- RunUMAP(HCC2_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC2_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC2_HPC,group.by = "group",label=T)




HCC3_HPC <- SCTransform(HCC3_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC3_HPC <- RunPCA(HCC3_HPC, verbose = FALSE)
HCC3_HPC <- RunUMAP(HCC3_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC3_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC3_HPC,group.by = "group",label=T)



HCC4_HPC <- SCTransform(HCC4_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC4_HPC <- RunPCA(HCC4_HPC, verbose = FALSE)
HCC4_HPC <- RunUMAP(HCC4_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC4_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC4_HPC,group.by = "group",label=T)




HCC5_HPC <- SCTransform(HCC5_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC5_HPC <- RunPCA(HCC5_HPC, verbose = FALSE)
HCC5_HPC <- RunUMAP(HCC5_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC5_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC5_HPC,group.by = "group",label=T)

HCC6_HPC <- SCTransform(HCC6_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC6_HPC <- RunPCA(HCC6_HPC, verbose = FALSE)
HCC6_HPC <- RunUMAP(HCC6_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC6_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC6_HPC,group.by = "group",label=T)




HCC7_HPC <- SCTransform(HCC7_HPC,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
HCC7_HPC <- RunPCA(HCC7_HPC, verbose = FALSE)
HCC7_HPC <- RunUMAP(HCC7_HPC, dims = 1:30, verbose = FALSE)
DimPlot(HCC7_HPC,group.by = "orig.ident",label=T)
DimPlot(HCC7_HPC,group.by = "group",label=T)





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


