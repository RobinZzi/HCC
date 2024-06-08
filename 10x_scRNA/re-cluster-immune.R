immune_merge
immune_merge <- SCTransform(immune_merge,method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
immune_merge <- RunPCA(immune_merge, verbose = FALSE)
immune_merge <- RunUMAP(immune_merge, dims = 1:30, verbose = FALSE)

DimPlot(immune_merge,group.by = "new.ident",label=T)
CairoPNG("FeaturePlot.png",xpd=TRUE,height = 800,width = 800)
FeaturePlot(immune_merge,features = "CD3D")
exha.markers <- FindMarkers(immune_merge,ident.1 = "CD8+ exhausted",ident.2 = "CD8+ cytotoxic")
