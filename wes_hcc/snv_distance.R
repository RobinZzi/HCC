setwd("~/projects/hcc/analysis/wes_hcc/hcc3")
rm(list=ls())
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
save.image("snv_distance.Rdata")


setwd("~/public_html/hcc_wes/hcc3/filter_gatk_result")
PT1 <-  fread("BCPT1_filter.vcf",skip = 264)

pt1_filt <- data.frame()
pt1_filt <- str_split_fixed(PT1$BCPT1,":",8)
colnames(pt1_filt) <- str_split_fixed(PT1[1,]$FORMAT,":",8)

pt1_filt <- as.data.frame(pt1_filt)
pt1_filt <- dplyr::select(pt1_filt,3,4)
pt1_filt <- lapply(pt1_filt, as.numeric)
pt1_filt <- as.data.frame(pt1_filt)
pt1_filt$sample <- 'PT1'

pt1_filt$pos <-  paste(PT1$`#CHROM`,PT1$POS,sep = "_")

ggplot(pt1_filt,aes(x=AF,y=DP))+geom_point()

pt1_filt


PT2 <-  fread("BCPT2_filter.vcf",skip = 264)

PT2_filt <- data.frame()
PT2_filt <- str_split_fixed(PT2$BCPT2,":",8)
colnames(PT2_filt) <- str_split_fixed(PT2[1,]$FORMAT,":",8)

PT2_filt <- as.data.frame(PT2_filt)
PT2_filt <- dplyr::select(PT2_filt,3,4)
PT2_filt <- lapply(PT2_filt, as.numeric)
PT2_filt <- as.data.frame(PT2_filt)

PT2_filt$sample <- 'PT2'
PT2_filt$pos <-  paste(PT2$`#CHROM`,PT2$POS,sep = "_")

ggplot(PT2_filt,aes(x=AF,y=DP))+geom_point()


PT3 <-  fread("BCPT3_filter.vcf",skip = 264)
PT3_filt <- data.frame()
PT3_filt <- str_split_fixed(PT3$BCPT3,":",8)
colnames(PT3_filt) <- str_split_fixed(PT3[1,]$FORMAT,":",8)

PT3_filt <- as.data.frame(PT3_filt)
PT3_filt <- dplyr::select(PT3_filt,3,4)
PT3_filt <- lapply(PT3_filt, as.numeric)
PT3_filt <- as.data.frame(PT3_filt)

PT3_filt$sample <- 'PT3'
PT3_filt$pos <-  paste(PT3$`#CHROM`,PT3$POS,sep = "_")

ggplot(PT3_filt,aes(x=AF,y=DP))+geom_point()




sum <- rbind(pt1_filt,PT2_filt,PT3_filt)

sum <- mutate(sum,mut_type = case_when(pos %in% hcc3_sum_pt1_pt2_table$pos ~ 'pt1_pt2',
                                       pos %in% hcc3_sum_pt1_pt3_table$pos ~ 'pt1_pt3',
                                       pos %in% hcc3_sum_pt2_pt3_table$pos ~ 'pt2_pt3',
                                       pos %in% hcc3_sum_pt1_pt2_pt3_table$pos ~ 'pt1_pt2_pt3') )


ggplot(sum,aes(x=AF,y=DP,color=mut_type))+geom_point()+ylim(0,500)



sum_subset <- subset(sum,subset= mut_type %in% c('pt1_pt2','pt1_pt3','pt2_pt3','pt1_pt2_pt3'))
sum_subset <- dplyr::select(sum_subset,4,3,5,1)
sum_subset_wide <- sum_subset %>%
  pivot_wider(names_from = sample, values_from = AF)



hcc2_sum <- dplyr::select(sum_subset_wide,3,4,5)
row.names(hcc2_sum) <- sum_subset_wide$pos
hcc2_sum[is.na(hcc2_sum)] <- 0



pheatmap(hcc2_sum,
         show_rownames = F,show_colnames = T,cluster_cols = T,cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,border_color = NA,
         treeheight_col = 10)


hcc2_sum_mtx <- t(hcc2_sum)
hcc2_sum_dist = dist(hcc2_sum_mtx, method = "euclidean")
hclust_hcc2_sum = hclust(hcc2_sum_dist, method = "complete")
plot(hclust_hcc2_sum)


hcc3_sum_pt1_pt2_table <- subset(hcc3_mut_table,subset = AAChange.refGene %in% row.names(hcc3_sum_pt1_pt2))
hcc3_sum_pt1_pt2_table$pos <- paste(hcc3_sum_pt1_pt2_table$Chromosome,hcc3_sum_pt1_pt2_table$Start_Position,sep = "_")
hcc3_sum_pt2_pt3_table <- subset(hcc3_mut_table,subset = AAChange.refGene %in% row.names(hcc3_sum_pt2_pt3))
hcc3_sum_pt2_pt3_table$pos <- paste(hcc3_sum_pt2_pt3_table$Chromosome,hcc3_sum_pt2_pt3_table$Start_Position,sep = "_")
hcc3_sum_pt1_pt3_table <- subset(hcc3_mut_table,subset = AAChange.refGene %in% row.names(hcc3_sum_pt1_pt3))
hcc3_sum_pt1_pt3_table$pos <- paste(hcc3_sum_pt1_pt3_table$Chromosome,hcc3_sum_pt1_pt3_table$Start_Position,sep = "_")
hcc3_sum_pt1_pt2_pt3_table <- subset(hcc3_mut_table,subset = AAChange.refGene %in% row.names(hcc3_sum_pt1_pt2_pt3))
hcc3_sum_pt1_pt2_pt3_table$pos <- paste(hcc3_sum_pt1_pt2_pt3_table$Chromosome,hcc3_sum_pt1_pt2_pt3_table$Start_Position,sep = "_")


hcc3_sum_subset <- rbind(hcc3_sum_pt1_pt2,hcc3_sum_pt1_pt3,hcc3_sum_pt2_pt3)

pheatmap(hcc3_sum_subset,
         show_rownames = F, show_colnames = T,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "complete",clustering_distance_cols = "canberra",
         color = colorRampPalette(c("#F3FBF2", "#AADCA9"))(20),
         treeheight_row = 0,
         treeheight_col = 12,border_color = NA,
         fontsize          = 12)
