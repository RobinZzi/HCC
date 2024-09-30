library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)

save.image("doublets.Rdata")

sweep.res.list <- paramSweep(bigseu, PCs = 1:20, sct = T)
# 对"scRNA_harmony"这个单细胞对象进行pN-pK参数扫描，以生成人工双细胞并计算每个细胞的pANN值
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
# 对汇总结果按照真实双胞胎比例（BCreal）进行升序排序，并显示排序后的数据框。这可以帮助我们找到双胞胎检测效果最好的参数组合
sweep.stats[order(sweep.stats$BCreal),]
# 根据汇总结果找到最优的pK参数。
bcmvn <- find.pK(sweep.stats) 
# 提取出全局最优的pK值，储存于"pK_bcmvn"
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
# 估计的同源双细胞（由两个或多个相同类型的细胞融合而成的假阳性细胞，它们通常比异源双细胞更难以检测和去除）的比例     
homotypic.prop <- modelHomotypic(bigseu$seurat_clusters) 
# 计算总的双细胞数量（假设双细胞形成率为 7.5%）
nExp_poi <- round(0.075 *nrow(bigseu@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 计算异源双细胞数量

# 使用确定好的参数鉴定doublets
bigseu_rmdl <- doubletFinder(bigseu, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                                  nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

# 可视化
DimPlot(bigseu_rmdl, group.by = "DF.classifications_0.25_26_6974", raster = FALSE)

# 将双细胞剔除后生成新的对象scRNA_harmony.singlet，便于我们后续分析
scRNA_harmony.singlet <- subset(scRNA_harmony, subset = DF.classifications_0.25_33_9717== "Singlet")