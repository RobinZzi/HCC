rm(list=ls())
setwd("~/projects/hcc/analysis/cell_re_rrbs/cnv_infer_norm")

save.image("CNV_norm.Rdata")

liver_bulk_1 <-fread("liver_bulk_1_10M.txt")
liver_bulk_1$x <- as.numeric(liver_bulk_1$x)
liver_bulk_2 <-fread("liver_bulk_2_10M.txt")
liver_bulk_2$x <- as.numeric(liver_bulk_2$x)

liver_bulk_1_sum <- sum(liver_bulk_1$x)
liver_bulk_2_sum <- sum(liver_bulk_2$x)

common_chr <- intersect(liver_bulk_1$Group.1,liver_bulk_2$Group.1)
liver_bulk_1_subset <- subset(liver_bulk_1,subset=Group.1 %in% common_chr)
liver_bulk_2_subset <- subset(liver_bulk_2,subset=Group.1 %in% common_chr)

liver_norm <- as.data.frame(cbind(liver_bulk_1_subset$x,liver_bulk_2_subset$x))
liver_norm$level <- rowMeans(liver_norm)
liver_norm$bin_name <- common_chr
liver_norm <- dplyr::select(liver_norm,4,3)




HepG2_scRRBS_1 <-fread("HepG2_scRRBS_1_10M.txt")
HepG2_scRRBS_1$x <- as.numeric(HepG2_scRRBS_1$x)
HepG2_scRRBS_2 <-fread("HepG2_scRRBS_2_10M.txt")
HepG2_scRRBS_2$x <- as.numeric(HepG2_scRRBS_2$x)

HepG2_scRRBS_1_sum <- sum(HepG2_scRRBS_1$x)
HepG2_scRRBS_2_sum <- sum(HepG2_scRRBS_2$x)
HepG2_scRRBS_1_co <- (HepG2_scRRBS_1_sum+HepG2_scRRBS_2_sum)/(2*HepG2_scRRBS_1_sum)
HepG2_scRRBS_2_co <- (HepG2_scRRBS_1_sum+HepG2_scRRBS_2_sum)/(2*HepG2_scRRBS_2_sum)


HepG2_scRRBS_1_aj <- HepG2_scRRBS_1
HepG2_scRRBS_1_aj$x <- HepG2_scRRBS_1_aj$x*HepG2_scRRBS_1_co
HepG2_scRRBS_1_aj <- dplyr::select(HepG2_scRRBS_1_aj,2,3)
colnames(HepG2_scRRBS_1_aj) <- c("bin_name","level")


HepG2_scRRBS_2_aj <- HepG2_scRRBS_2
HepG2_scRRBS_2_aj$x <- HepG2_scRRBS_2_aj$x*HepG2_scRRBS_2_co
HepG2_scRRBS_2_aj <- dplyr::select(HepG2_scRRBS_2_aj,2,3)
colnames(HepG2_scRRBS_2_aj) <- c("bin_name","level")



HepG2_bulk_RRBS <-fread("HepG2_bulk_RRBS_10M.txt")

HepG2_bulk_RRBS_sub <- subset(HepG2_bulk_RRBS,subset=Group.1 %in% common_chr)

HepG2_bulk_RRBS_sub <- dplyr::select(HepG2_bulk_RRBS_sub,2,3)
colnames(HepG2_bulk_RRBS_sub) <- c("bin_name","level")





HepG2_scTrio_1 <-fread("HepG2_scTrio_1_10M.txt")
HepG2_scTrio_1$x <- as.numeric(HepG2_scTrio_1$x)
HepG2_scTrio_2 <-fread("HepG2_scTrio_2_10M.txt")
HepG2_scTrio_2$x <- as.numeric(HepG2_scTrio_2$x)
HepG2_scTrio_3 <-fread("HepG2_scTrio_3_10M.txt")
HepG2_scTrio_3$x <- as.numeric(HepG2_scTrio_3$x)
HepG2_scTrio_4 <-fread("HepG2_scTrio_4_10M.txt")
HepG2_scTrio_4$x <- as.numeric(HepG2_scTrio_4$x)
HepG2_scTrio_5 <-fread("HepG2_scTrio_5_10M.txt")
HepG2_scTrio_5$x <- as.numeric(HepG2_scTrio_5$x)
HepG2_scTrio_6 <-fread("HepG2_scTrio_6_10M.txt")
HepG2_scTrio_6$x <- as.numeric(HepG2_scTrio_6$x)


HepG2_scTrio_1_sum <- sum(HepG2_scTrio_1$x)
HepG2_scTrio_2_sum <- sum(HepG2_scTrio_2$x)
HepG2_scTrio_3_sum <- sum(HepG2_scTrio_3$x)
HepG2_scTrio_4_sum <- sum(HepG2_scTrio_4$x)
HepG2_scTrio_5_sum <- sum(HepG2_scTrio_5$x)
HepG2_scTrio_6_sum <- sum(HepG2_scTrio_6$x)

HepG2_scTrio_1_co <- (HepG2_scTrio_1_sum+HepG2_scTrio_2_sum+HepG2_scTrio_3_sum+
                      HepG2_scTrio_4_sum+HepG2_scTrio_5_sum+HepG2_scTrio_6_sum)/(6*HepG2_scTrio_1_sum)
HepG2_scTrio_2_co <- (HepG2_scTrio_1_sum+HepG2_scTrio_2_sum+HepG2_scTrio_3_sum+
                      HepG2_scTrio_4_sum+HepG2_scTrio_5_sum+HepG2_scTrio_6_sum)/(6*HepG2_scTrio_2_sum)
HepG2_scTrio_3_co <- (HepG2_scTrio_1_sum+HepG2_scTrio_2_sum+HepG2_scTrio_3_sum+
                        HepG2_scTrio_4_sum+HepG2_scTrio_5_sum+HepG2_scTrio_6_sum)/(6*HepG2_scTrio_3_sum)
HepG2_scTrio_4_co <- (HepG2_scTrio_1_sum+HepG2_scTrio_2_sum+HepG2_scTrio_3_sum+
                        HepG2_scTrio_4_sum+HepG2_scTrio_5_sum+HepG2_scTrio_6_sum)/(6*HepG2_scTrio_4_sum)
HepG2_scTrio_5_co <- (HepG2_scTrio_1_sum+HepG2_scTrio_2_sum+HepG2_scTrio_3_sum+
                        HepG2_scTrio_4_sum+HepG2_scTrio_5_sum+HepG2_scTrio_6_sum)/(6*HepG2_scTrio_5_sum)
HepG2_scTrio_6_co <- (HepG2_scTrio_1_sum+HepG2_scTrio_2_sum+HepG2_scTrio_3_sum+
                        HepG2_scTrio_4_sum+HepG2_scTrio_5_sum+HepG2_scTrio_6_sum)/(6*HepG2_scTrio_6_sum)


HepG2_scTrio_1_aj <- HepG2_scTrio_1
HepG2_scTrio_1_aj$x <- HepG2_scTrio_1_aj$x*HepG2_scTrio_1_co
HepG2_scTrio_1_aj <- dplyr::select(HepG2_scTrio_1_aj,2,3)
colnames(HepG2_scTrio_1_aj) <- c("bin_name","level")

HepG2_scTrio_2_aj <- HepG2_scTrio_2
HepG2_scTrio_2_aj$x <- HepG2_scTrio_2_aj$x*HepG2_scTrio_2_co
HepG2_scTrio_2_aj <- dplyr::select(HepG2_scTrio_2_aj,2,3)
colnames(HepG2_scTrio_2_aj) <- c("bin_name","level")

HepG2_scTrio_3_aj <- HepG2_scTrio_3
HepG2_scTrio_3_aj$x <- HepG2_scTrio_3_aj$x*HepG2_scTrio_3_co
HepG2_scTrio_3_aj <- dplyr::select(HepG2_scTrio_3_aj,2,3)
colnames(HepG2_scTrio_3_aj) <- c("bin_name","level")

HepG2_scTrio_4_aj <- HepG2_scTrio_4
HepG2_scTrio_4_aj$x <- HepG2_scTrio_4_aj$x*HepG2_scTrio_4_co
HepG2_scTrio_4_aj <- dplyr::select(HepG2_scTrio_4_aj,2,3)
colnames(HepG2_scTrio_4_aj) <- c("bin_name","level")

HepG2_scTrio_5_aj <- HepG2_scTrio_5
HepG2_scTrio_5_aj$x <- HepG2_scTrio_5_aj$x*HepG2_scTrio_5_co
HepG2_scTrio_5_aj <- dplyr::select(HepG2_scTrio_5_aj,2,3)
colnames(HepG2_scTrio_5_aj) <- c("bin_name","level")

HepG2_scTrio_6_aj <- HepG2_scTrio_6
HepG2_scTrio_6_aj$x <- HepG2_scTrio_6_aj$x*HepG2_scTrio_6_co
HepG2_scTrio_6_aj <- dplyr::select(HepG2_scTrio_6_aj,2,3)
colnames(HepG2_scTrio_6_aj) <- c("bin_name","level")



all_common_chr <- Reduce(intersect,list(liver_norm$bin_name,HepG2_scRRBS_1_aj$bin_name,HepG2_scRRBS_2_aj$bin_name,
                        HepG2_bulk_RRBS_sub$bin_name,
                        HepG2_scTrio_1_aj$bin_name,HepG2_scTrio_2_aj$bin_name,HepG2_scTrio_3_aj$bin_name,
                        HepG2_scTrio_4_aj$bin_name,HepG2_scTrio_5_aj$bin_name,HepG2_scTrio_6_aj$bin_name))



liver_norm <- subset(liver_norm,subset=bin_name %in% all_common_chr)

HepG2_bulk_RRBS_sub <- subset(HepG2_bulk_RRBS_sub,subset=bin_name %in% all_common_chr)

HepG2_scRRBS_1_aj_sub <- subset(HepG2_scRRBS_1_aj,subset=bin_name %in% all_common_chr)
HepG2_scRRBS_2_aj_sub <- subset(HepG2_scRRBS_2_aj,subset=bin_name %in% all_common_chr)

HepG2_scTrio_1_aj_sub <- subset(HepG2_scTrio_1_aj,subset=bin_name %in% all_common_chr)
HepG2_scTrio_2_aj_sub <- subset(HepG2_scTrio_2_aj,subset=bin_name %in% all_common_chr)
HepG2_scTrio_3_aj_sub <- subset(HepG2_scTrio_3_aj,subset=bin_name %in% all_common_chr)
HepG2_scTrio_4_aj_sub <- subset(HepG2_scTrio_4_aj,subset=bin_name %in% all_common_chr)
HepG2_scTrio_5_aj_sub <- subset(HepG2_scTrio_5_aj,subset=bin_name %in% all_common_chr)
HepG2_scTrio_6_aj_sub <- subset(HepG2_scTrio_6_aj,subset=bin_name %in% all_common_chr)



sum_mtx <- as.data.frame(cbind(HepG2_bulk_RRBS_sub$level,HepG2_scRRBS_1_aj_sub$level,HepG2_scRRBS_2_aj_sub$level,
                               HepG2_scTrio_1_aj_sub$level,HepG2_scTrio_2_aj_sub$level,HepG2_scTrio_3_aj_sub$level,
                               HepG2_scTrio_4_aj_sub$level,HepG2_scTrio_5_aj_sub$level,HepG2_scTrio_6_aj_sub$level))

sum_mtx_aj <- sweep(sum_mtx,1, liver_norm$level, `/`)  

colnames(sum_mtx_aj) <- c("bulk_RRBS","HepG2_scRRBS_1","HepG2_scRRBS_2",
                          "HepG2_scTrio_1","HepG2_scTrio_2","HepG2_scTrio_3",
                          "HepG2_scTrio_4","HepG2_scTrio_5","HepG2_scTrio_6")

row.names(sum_mtx_aj) <- all_common_chr

sum_mtx_aj_t <- as.data.frame(t(sum_mtx_aj))
col_anno <- as.data.frame(str_split_fixed(all_common_chr,"_",2))
col_anno[,3:4] <- str_split_fixed(col_anno$V1,"r",2)
col_anno$V2 <- as.numeric(col_anno$V2)
col_anno$V4 <- as.numeric(col_anno$V4)
col_anno <- mutate(col_anno,id = 100*V4+V2)
row.names(col_anno) <- all_common_chr
col_anno <- col_anno[order(col_anno$id),]
row_names_factor <- factor(rownames(col_anno)) 
col_anno$bin_name <- rownames(col_anno)
levels(row_names_factor) <- col_anno$bin_name

colnames(sum_mtx_aj_t) <- factor(colnames(sum_mtx_aj_t),levels = levels(row_names_factor))
col_anno_2 <- as.data.frame(dplyr::select(col_anno,1))
colnames(col_anno_2) <- c("chr")

sum_mtx_aj_t_sort <- dplyr::select(sum_mtx_aj_t,levels(row_names_factor))

col_anno_2$chr <- factor(col_anno_2$chr,levels = c("chr1","chr2","chr3","chr4","chr5",
                                                   "chr6","chr7","chr8","chr9","chr10",
                                                   "chr11","chr12","chr13","chr14","chr15",
                                                   "chr16","chr17","chr18","chr19","chr20",
                                                   "chr21","chr22","chrX","chrY"))



ann_colors =list(
  chr=c('chr1'='#E5D2DD','chr2'='#53A85F','chr3'= '#F1BB72','chr4'='#F3B1A0','chr5'='#D6E7A3',
        'chr6'='#57C3F3','chr7'='#476D87','chr8'='#E95C59','chr9'= '#E59CC4','chr10'= '#AB3282', 
        'chr11'='#23452F', 'chr12'='#BD956A', 'chr13'='#8C549C', 'chr14'='#585658', 'chr15'='#9FA3A8',
        'chr16'='#E0D4CA', 'chr17'='#5F3D69','chr18'= '#C5DEBA','chr19'= '#58A4C3','chr20'='#E4C755', 
        'chr21'='#F7F398', 'chr22'='#AA9A59','chrX'= '#E63863', 'chrY'='#E39A35')
)

order_row <- c("HepG2_scTrio_1","HepG2_scTrio_2","HepG2_scTrio_3",
               "HepG2_scTrio_4","HepG2_scTrio_5","HepG2_scTrio_6",
               "HepG2_scRRBS_1","HepG2_scRRBS_2","bulk_RRBS")

pheatmap(sum_mtx_aj_t_sort[order_row,],
         scale = 'row',cluster_rows = F,cluster_cols = F,angle_col = 315,breaks = seq(-2, 2, length.out = 100),
         show_colnames = F,annotation_col = col_anno_2,annotation_colors = ann_colors,
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         cellwidth = 1,cellheight = 10,ylab = "left")


col_anno_2_chr1418 <- subset(col_anno_2 ,subset = chr %in% c("chr14","chr15","chr16","chr17","chr18"))
sum_mtx_aj_t_sort_1418 <- dplyr::select(sum_mtx_aj_t_sort[order_row,],row.names(col_anno_2_chr1418))
pheatmap(sum_mtx_aj_t_sort_1418,
         scale = 'row',cluster_rows = F,cluster_cols = F,angle_col = 315,breaks = seq(-2, 2, length.out = 100),
         show_colnames = F,annotation_col = col_anno_2_chr1418,annotation_colors = ann_colors,
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         cellwidth = 10,cellheight = 12,ylab = "left")
