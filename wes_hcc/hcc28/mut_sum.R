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
library(rjags)

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
setwd("~/projects/hcc/analysis/wes_hcc")
load("hcc28.Rdata")
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
write.table(hcc3_mut_table,"hcc3_mut_table.txt")





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
oncoplot(maf = hcc3_laml,
         top = 20,
         fontSize = 0.5,
         showTumorSampleBarcodes = T)
dev.off()
maftools::lollipopPlot(hcc28_laml,
                       gene = 'TTN',
                       AACol = 'AAChange.refGene',
                       labelPos = 'all')
maftools::lollipopPlot(hcc3_laml,
                       gene = 'TTN',
                       AACol = 'AAChange.refGene',
                       labelPos = 'all')
maftools::lollipopPlot(hcc28_laml,
                       gene = 'NNT',
                       AACol = 'AAChange.refGene',
                       labelPos = 'all')
maftools::lollipopPlot(hcc29_laml,
                       gene = 'NNT',
                       AACol = 'AAChange.refGene',
                       labelPos = 'all')

hcc3_mut_table <- as.data.frame(hcc3_laml@data)
hcc3_mut_table$Tumor_Sample_Barcode
hcc3_mut_table$Tumor_Sample_Barcode[is.na(hcc3_mut_table$Tumor_Sample_Barcode)] <- 'PT4'
for (i in 1:nrow(hcc3_mut_table)) {
  
  if(hcc3_mut_table[i,]$Tumor_Sample_Barcode %in% c("PT1","PT2","PT3")){
  
  }else{
    hcc3_mut_table[i,]$Tumor_Sample_Barcode <- 'PT4'
  }
}


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
sum_hcc29 <- intersect(unique(hcc29_pt1_mut_table$Hugo_Symbol),unique(hcc29_pt4_mut_table$Hugo_Symbol))


sum3_28 <- intersect(sum_hcc3,sum_hcc28)
sum29_28 <- intersect(sum_hcc29,sum_hcc28)
sum29_3 <- intersect(sum_hcc29,sum_hcc3)



hcc28_common_mut_go <- enrichGO(gene  = sum_hcc28,
                                           OrgDb      = org.Hs.eg.db,
                                           keyType    = 'SYMBOL',
                                           ont        = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 0.05,
                                           qvalueCutoff = 0.05)
hcc28_common_mut_go <- as.data.frame(hcc28_common_mut_go@result)
hcc28_common_mut_go [,"logp"] <- -log10(hcc28_common_mut_go$pvalue)
ggplot(data = hcc28_common_mut_go[1:20,])+
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
       title="Common mut in hcc28")

hcc29_common_mut_go <- enrichGO(gene  = sum_hcc29,
                                OrgDb      = org.Hs.eg.db,
                                keyType    = 'SYMBOL',
                                ont        = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
hcc29_common_mut_go <- as.data.frame(hcc29_common_mut_go@result)
hcc29_common_mut_go [,"logp"] <- -log10(hcc29_common_mut_go$pvalue)
ggplot(data = hcc29_common_mut_go[1:20,])+
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
       title="Common mut in hcc29")




hcc3_common_mut_go <- enrichGO(gene  = sum_hcc3,
                                OrgDb      = org.Hs.eg.db,
                                keyType    = 'SYMBOL',
                                ont        = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
hcc3_common_mut_go <- as.data.frame(hcc3_common_mut_go@result)
hcc3_common_mut_go [,"logp"] <- -log10(hcc3_common_mut_go$pvalue)
ggplot(data = hcc3_common_mut_go[1:20,])+
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
       title="Common mut in hcc29")


hcc28_info <- select(hcc28_mut_table,Tumor_Sample_Barcode,AAChange.refGene)

hcc28_info_PT1 <- subset(hcc28_info,subset = hcc28_info$Tumor_Sample_Barcode =='PT1')
hcc28_info_PT1 <- hcc28_info_PT1[1:165,]
row.names(hcc28_info_PT1) <- hcc28_info_PT1$AAChange.refGene
colnames(hcc28_info_PT1) <- c("PT1","pos")
hcc28_info_PT1$PT1 <- 1

hcc28_info_PT2 <- subset(hcc28_info,subset = hcc28_info$Tumor_Sample_Barcode =='PT2')
hcc28_info_PT2 <- hcc28_info_PT2[1:130,]
row.names(hcc28_info_PT2) <- hcc28_info_PT2$AAChange.refGene
colnames(hcc28_info_PT2) <- c("PT2","pos")
hcc28_info_PT2$PT2 <- 1

hcc28_info_PT4 <- subset(hcc28_info,subset = hcc28_info$Tumor_Sample_Barcode =='PT4')
hcc28_info_PT4 <- hcc28_info_PT4[1:141,]
row.names(hcc28_info_PT4) <- hcc28_info_PT4$AAChange.refGene
colnames(hcc28_info_PT4) <- c("PT4","pos")
hcc28_info_PT4$PT4 <- 1




hcc28_sum <- merge(hcc28_info_PT1,hcc28_info_PT2,all=T)
row.names(hcc28_sum) <- hcc28_sum$pos

hcc28_sum <- merge(hcc28_sum,hcc28_info_PT4,all=T)
hcc28_sum[is.na(hcc28_sum)] <- 0
row.names(hcc28_sum) <- hcc28_sum$pos
hcc28_sum <- hcc28_sum[,-1]
pheatmap(hcc28_sum,
         show_rownames = F, show_colnames = T,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",
         color = colorRampPalette(c("#fff5f6", "#fc929f"))(20),
         treeheight_row = 0,
         treeheight_col = 20,
         fontsize          = 12)


hcc28_sum_hclust <- as.data.frame(t(hcc28_sum))
hcc28_mut_dist = dist(hcc28_sum_hclust, method = "euclidean")
hclust_hcc28 = hclust(hcc28_mut_dist, method = "average")
plot(hclust_hcc28)













hcc29_test <- unique(hcc29_mut_table$AAChange.refGene)
hcc29_info <- select(hcc29_mut_table,Tumor_Sample_Barcode,AAChange.refGene)

hcc29_info_PT1 <- subset(hcc29_info,subset = hcc29_info$Tumor_Sample_Barcode =='PT1')
row.names(hcc29_info_PT1) <- hcc29_info_PT1$AAChange.refGene
colnames(hcc29_info_PT1) <- c("PT1","pos")
hcc29_info_PT1$PT1 <- 1

hcc29_info_PT3 <- subset(hcc29_info,subset = hcc29_info$Tumor_Sample_Barcode =='PT3')
row.names(hcc29_info_PT3) <- hcc29_info_PT3$AAChange.refGene
colnames(hcc29_info_PT3) <- c("PT3","pos")
hcc29_info_PT3$PT3 <- 1

hcc29_info_PT4 <- subset(hcc29_info,subset = hcc29_info$Tumor_Sample_Barcode =='PT4')
row.names(hcc29_info_PT4) <- hcc29_info_PT4$AAChange.refGene
colnames(hcc29_info_PT4) <- c("PT4","pos")
hcc29_info_PT4$PT4 <- 1




hcc29_sum <- merge(hcc29_info_PT1,hcc29_info_PT3,all=T)
row.names(hcc29_sum) <- hcc29_sum$pos

hcc29_sum <- merge(hcc29_sum,hcc29_info_PT4,all=T)
hcc29_sum[is.na(hcc29_sum)] <- 0
row.names(hcc29_sum) <- hcc29_sum$pos
hcc29_sum <- hcc29_sum[,-1]
pheatmap(hcc29_sum,
         show_rownames = F, show_colnames = T,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",
         color = colorRampPalette(c("#fff5f6", "#fc929f"))(20),
         treeheight_row = 0,
         treeheight_col = 20,border_color = NA,
         fontsize          = 12)


hcc29_sum_hclust <- as.data.frame(t(hcc29_sum))
hcc29_mut_dist = dist(hcc29_sum_hclust, method = "euclidean")
hclust_hcc29 = hclust(hcc29_mut_dist, method = "average")
plot(hclust_hcc29)








hcc3_mut_table <- mut_table
hcc3_info <- select(hcc3_mut_table,Tumor_Sample_Barcode,AAChange.refGene)

hcc3_info_PT1 <- subset(hcc3_info,subset = hcc3_info$Tumor_Sample_Barcode =='BCPT1')
row.names(hcc3_info_PT1) <- hcc3_info_PT1$AAChange.refGene
colnames(hcc3_info_PT1) <- c("PT1","pos")
hcc3_info_PT1$PT1 <- 1

hcc3_info_PT2 <- subset(hcc3_info,subset = hcc3_info$Tumor_Sample_Barcode =='BCPT2')
row.names(hcc3_info_PT2) <- hcc3_info_PT3$AAChange.refGene
colnames(hcc3_info_PT2) <- c("PT2","pos")
hcc3_info_PT2$PT2 <- 1

hcc3_info_PT3 <- subset(hcc3_info,subset = hcc3_info$Tumor_Sample_Barcode =='BCPT3')
row.names(hcc3_info_PT3) <- hcc3_info_PT3$AAChange.refGene
colnames(hcc3_info_PT3) <- c("PT3","pos")
hcc3_info_PT3$PT3 <- 1

hcc3_info_PT4 <- subset(hcc3_info,subset = hcc3_info$Tumor_Sample_Barcode =='BCPT4')
row.names(hcc3_info_PT4) <- hcc3_info_PT4$AAChange.refGene
colnames(hcc3_info_PT4) <- c("PT4","pos")
hcc3_info_PT4$PT4 <- 1




hcc3_sum <- merge(hcc3_info_PT1,hcc3_info_PT3,all=T)
row.names(hcc3_sum) <- hcc3_sum$pos

hcc3_sum <- merge(hcc3_sum,hcc3_info_PT4,all=T)
hcc3_sum <- merge(hcc3_sum,hcc3_info_PT2,all=T)


hcc3_sum[is.na(hcc3_sum)] <- 0
hcc3_sum <- hcc3_sum[-1,]
row.names(hcc3_sum) <- hcc3_sum$pos
hcc3_sum <- hcc3_sum[,-1]
pheatmap(hcc3_sum,
         show_rownames = F, show_colnames = T,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "complete",clustering_distance_cols = "canberra",
         color = colorRampPalette(c("#fff5f6", "#fc929f"))(20),
         treeheight_row = 0,
         treeheight_col = 20,border_color = NA,
         fontsize          = 12)


hcc3_sum_hclust <- as.data.frame(t(hcc3_sum))
hcc3_mut_dist = dist(hcc3_sum_hclust, method = "canberra")
hclust_hcc3 = hclust(hcc3_mut_dist, method = "average")
plot(hclust_hcc3)


