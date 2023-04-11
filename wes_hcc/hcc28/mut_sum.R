library('maftools')
library(BiocManager)
library(maftools)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)
library(bedtoolsr)
library(org.Hs.eg.db)
library(ggplot2)
library(ggVennDiagram)
library(venn)
library(VennDiagram)


options(stringsAsFactors = F)

rm(list = ls())
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc29/gatk_result")
save.image("mut_sum.Rdata")
hcc29.annovar.maf <- annovarToMaf(annovar = "hcc29.maf",
                                refBuild = 'hg38',
                                tsbCol = 'Tumor_Sample_Barcode1',
                                table = 'refGene',
                                MAFobj = T)
hcc29_laml <- hcc29.annovar.maf



unique(hcc29_laml@data$Tumor_Sample_Barcode)
getSampleSummary(hcc29_laml)
getGeneSummary(hcc29_laml)
getFields(hcc29_laml)




png(paste0('plotmafSummary_', project, '.png'),res = 150,width = 1080,height = 1080)
plotmafSummary(maf = hcc29_laml,
               rmOutlier = TRUE,
               showBarcodes = T,
               textSize = 0.4,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE)
dev.off()
png(paste0('oncoplot_top30_', project, '.png'),res = 150,width = 1080,height = 1080)
oncoplot(maf = hcc29_laml,
         top = 20,
         fontSize = 0.5,
         showTumorSampleBarcodes = T)
dev.off()
maftools::lollipopPlot(hcc28_laml,
                       gene = 'TP53',
                       AACol = 'AAChange.refGene',
                       labelPos = 'all')


hcc29_mut_table <- as.data.frame(hcc29_laml@data)





setwd("~/projects/hcc/scripts/wes_hcc/hcc28")
library('maftools')
library(BiocManager)
library(maftools)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)
library(bedtoolsr)
library(org.Hs.eg.db)
library(ggplot2)
library(ggVennDiagram)
library(venn)
options(stringsAsFactors = F)

save.image("hcc28.Rdata")
rm(list = ls())
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc28/gatk_result")
hcc28.annovar.maf <- annovarToMaf(annovar = "hcc28.maf",
                                  refBuild = 'hg38',
                                  tsbCol = 'Tumor_Sample_Barcode1',
                                  table = 'refGene',
                                  MAFobj = T)
hcc28_laml <- hcc28.annovar.maf



unique(hcc28_laml@data$Tumor_Sample_Barcode)
getSampleSummary(hcc28_laml)
getGeneSummary(hcc28_laml)
getFields(hcc28_laml)




png(paste0('plotmafSummary_', project, '.png'),res = 150,width = 1080,height = 1080)
plotmafSummary(maf = hcc28_laml,
               rmOutlier = TRUE,
               showBarcodes = T,
               textSize = 0.4,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE)
dev.off()
png(paste0('oncoplot_top30_', project, '.png'),res = 150,width = 1080,height = 1080)
oncoplot(maf = hcc28_laml,
         top = 20,
         fontSize = 0.5,
         showTumorSampleBarcodes = T)
dev.off()
maftools::lollipopPlot(hcc28_laml,
                       gene = 'CTNNB1',
                       AACol = 'AAChange.refGene',
                       labelPos = 'all')


hcc28_mut_table <- as.data.frame(hcc28_laml@data)




write.table(hcc28_mut_table,"hcc28_mut_table.txt")
write.table(hcc29_mut_table,"hcc29_mut_table.txt")






hcc28_pt1_mut_table <- subset(hcc28_mut_table,subset= hcc28_mut_table$Tumor_Sample_Barcode == "PT1")
hcc28_pt2_mut_table <- subset(hcc28_mut_table,subset= hcc28_mut_table$Tumor_Sample_Barcode == "PT2")
hcc28_pt4_mut_table <- subset(hcc28_mut_table,subset= hcc28_mut_table$Tumor_Sample_Barcode == "PT4")

hcc28_pt1_unique_mut <- hcc28_pt1_mut_table$AAChange.refGene #86
hcc28_pt2_unique_mut <- hcc28_pt2_mut_table$AAChange.refGene #86
hcc28_pt4_unique_mut <- hcc28_pt4_mut_table$AAChange.refGene #56



hcc29_pt1_mut_table <- subset(hcc29_mut_table,subset= hcc29_mut_table$Tumor_Sample_Barcode == "PT1")
hcc29_pt3_mut_table <- subset(hcc29_mut_table,subset= hcc29_mut_table$Tumor_Sample_Barcode == "PT3")
hcc29_pt4_mut_table <- subset(hcc29_mut_table,subset= hcc29_mut_table$Tumor_Sample_Barcode == "PT4")

hcc29_pt1_unique_mut <- hcc29_pt1_mut_table$AAChange.refGene#86
hcc29_pt3_unique_mut <- hcc29_pt3_mut_table$AAChange.refGene #86
hcc29_pt4_unique_mut <- hcc29_pt4_mut_table$AAChange.refGene #56

hcc28_Venn <- list(hcc28_pt1_unique_mut,hcc28_pt2_unique_mut,hcc28_pt4_unique_mut)
hcc29_Venn <- list(hcc29_pt1_unique_mut,hcc29_pt3_unique_mut,hcc29_pt4_unique_mut)
venn(hcc28_Venn,zcolor=c('red','yellow','blue','green'),ilcs=0.)

hcc28_Venn <- list(pt1 = hcc28_pt1_unique_mut, pt2 = hcc28_pt2_unique_mut, pt4 = hcc28_pt4_unique_mut)
hcc29_Venn <- list(pt1 = hcc29_pt1_unique_mut, pt3 = hcc29_pt3_unique_mut, pt4 = hcc29_pt4_unique_mut)


venn.diagram(hcc28_Venn, filename = 'hcc28.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

venn.diagram(hcc29_Venn, filename = 'hcc29.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')















setwd("~/projects/hcc/scripts/wes_hcc/hcc3")
library('maftools')
library(BiocManager)
library(maftools)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)
library(bedtoolsr)
library(org.Hs.eg.db)
library(ggplot2)
library(ggVennDiagram)
library(venn)
options(stringsAsFactors = F)

save.image("hcc3.Rdata")
rm(list = ls())
setwd("~/projects/hcc/analysis/wes_hcc/hcc3")
hcc3.annovar.maf <- annovarToMaf(annovar = "hcc3.maf",
                                  refBuild = 'hg38',
                                  tsbCol = 'Tumor_Sample_Barcode1',
                                  table = 'refGene',
                                  MAFobj = T)
hcc3_laml <- hcc3.annovar.maf



unique(hcc3_laml@data$Tumor_Sample_Barcode)
getSampleSummary(hcc3_laml)
getGeneSummary(hcc3_laml)
getFields(hcc3_laml)




png(paste0('plotmafSummary_', project, '.png'),res = 150,width = 1080,height = 1080)
plotmafSummary(maf = hcc3_laml,
               rmOutlier = TRUE,
               showBarcodes = T,
               textSize = 0.4,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE)
dev.off()
png(paste0('oncoplot_top30_', project, '.png'),res = 150,width = 1080,height = 1080)
oncoplot(maf = hcc28_laml,
         top = 20,
         fontSize = 0.5,
         showTumorSampleBarcodes = T)
dev.off()
maftools::lollipopPlot(hcc28_laml,
                       gene = 'CTNNB1',
                       AACol = 'AAChange.refGene',
                       labelPos = 'all')


hcc3_mut_table <- as.data.frame(hcc3_laml@data)




write.table(hcc3_mut_table,"hcc3_mut_table.txt")






hcc3_pt1_mut_table <- subset(hcc3_mut_table,subset= hcc3_mut_table$Tumor_Sample_Barcode == "PT1")
hcc3_pt2_mut_table <- subset(hcc3_mut_table,subset= hcc3_mut_table$Tumor_Sample_Barcode == "PT2")
hcc3_pt3_mut_table <- subset(hcc3_mut_table,subset= hcc3_mut_table$Tumor_Sample_Barcode == "PT3")
hcc3_pt4_mut_table <- subset(hcc3_mut_table,subset= hcc3_mut_table$Tumor_Sample_Barcode == "PT4")


hcc3_pt1_unique_mut <- hcc3_pt1_mut_table$AAChange.refGene #86
hcc3_pt2_unique_mut <- hcc3_pt2_mut_table$AAChange.refGene #86
hcc3_pt3_unique_mut <- hcc3_pt3_mut_table$AAChange.refGene #86
hcc3_pt4_unique_mut <- hcc3_pt4_mut_table$AAChange.refGene #56



hcc3_Venn <- list(hcc3_pt1_unique_mut,hcc3_pt2_unique_mut,hcc3_pt3_unique_mut)

hcc28_Venn <- list(hcc28_pt1_unique_mut,hcc28_pt2_unique_mut,hcc28_pt4_unique_mut)
hcc29_Venn <- list(hcc29_pt1_unique_mut,hcc29_pt3_unique_mut,hcc29_pt4_unique_mut)
venn(hcc28_Venn,zcolor=c('red','yellow','blue','green'),ilcs=0.)



hcc3_Venn <- list(pt1 = hcc3_pt1_unique_mut, pt2 = hcc3_pt2_unique_mut, pt3 = hcc3_pt3_unique_mut)
hcc28_Venn <- list(pt1 = hcc28_pt1_unique_mut, pt2 = hcc28_pt2_unique_mut, pt4 = hcc28_pt4_unique_mut)
hcc29_Venn <- list(pt1 = hcc29_pt1_unique_mut, pt3 = hcc29_pt3_unique_mut, pt4 = hcc29_pt4_unique_mut)

sum_Venn <- list(hcc3_pt1 = hcc3_pt1_unique_mut, hcc3_pt2 = hcc3_pt2_unique_mut, hcc3_pt3 = hcc3_pt3_unique_mut,hcc28_pt1 = hcc28_pt1_unique_mut, hcc28_pt2 = hcc28_pt2_unique_mut, hcc28_pt4 = hcc28_pt4_unique_mut,hcc29_pt1 = hcc29_pt1_unique_mut, hcc29_pt3 = hcc29_pt3_unique_mut, hcc29_pt4 = hcc29_pt4_unique_mut)

venn.diagram(hcc28_Venn, filename = 'hcc28.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

venn.diagram(hcc29_Venn, filename = 'hcc29.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')


venn.diagram(hcc3_Venn, filename = 'hcc3.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')



sum_hcc3 <- intersect(unique(hcc3_pt1_mut_table$Hugo_Symbol),unique(hcc3_pt2_mut_table$Hugo_Symbol))
sum_hcc28 <- intersect(unique(hcc28_pt1_mut_table$Hugo_Symbol),unique(hcc28_pt2_mut_table$Hugo_Symbol))
sum_hcc29 <- intersect(unique(hcc28_pt1_mut_table$Hugo_Symbol),unique(hcc29_pt3_mut_table$Hugo_Symbol))


sum3_28 <- intersect(sum_hcc3,sum_hcc28)
sum29_28 <- intersect(sum_hcc29,sum_hcc28)
sum29_3 <- intersect(sum_hcc29,sum_hcc28)
