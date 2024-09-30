setwd("~/projects/hcc/analysis/trio_seq/bismark/cnvkit_result/scMerge/heatmap")
rm()

library(data.table)
library(dplyr)  
library(readr)  
library(tidyverse)
library(pheatmap)  
save.image("cnvkit_distance.Rdata")

# 步骤 1: 读取所有 cnr 文件  
setwd("~/projects/hcc/analysis/trio_seq/bismark/cnvkit_result/scMerge/sum")
path_to_cnr_files <- getwd() 
file_list <- list.files(path_to_cnr_files, pattern="*.cns", full.names=TRUE)  

# 初始化一个空的数据框来存储合并的数据  
combined_df <- data.frame()  

# 读取每个文件并合并数据  
for (file in file_list) {  
  sample_name <- tools::file_path_sans_ext(basename(file))  
  df <- read_tsv(file, col_types = cols())  
  if(is.null(combined_df[["chromosome"]])) {  # 只在第一次迭代时加入位置列  
    combined_df <- df %>%  
      select(chromosome, start, end, log2) %>%  
      rename(!!sample_name := log2)  
  } else {  
    combined_df <- combined_df %>%  
      mutate(!!sample_name := df$log2)  
  }  
}  

pq_info <- c(rep(c('p','q'),12),"p",rep(c('p','q'),4),"p",rep(c('p','q'),2),"p",rep(c('p','q'),2))
combined_df$pq <- pq_info
row.names(combined_df) <- paste(combined_df$chromosome,combined_df$pq,sep = "_")

# 步骤 2: 移除不需要的列，仅保留 log2 ratios  
combined_df_mtx <- combined_df %>% select(-chromosome,-start,-end,-pq)  
row.names(combined_df_mtx) <- row.names(combined_df)


combined_df_mtx_fix <- rbind(combined_df_mtx,rep(0,55),rep(0,55),rep(0,55))


row.names(combined_df_mtx_fix) <- c(row.names(combined_df),"chr13_q","chr18_q","chr21_q")

agg_mtx <- t(combined_df_mtx_fix)
pt1_mtx <- agg_mtx[1:34,]
pt2_mtx <- agg_mtx[35:50,]
pt3_mtx <- agg_mtx[51:55,]



agged_mtx <- rbind(colMeans(pt1_mtx),colMeans(pt2_mtx),colMeans(pt3_mtx))
row.names(agged_mtx) <- c('pt1','pt2','pt3')
scaled_agged_mtx <- scale(agged_mtx)

# 步骤 3: 填补缺失值并标准化数据  
filled_df <- combined_df_mtx_fix %>% replace(is.na(.), 0)  
standardized_df <- scale(t(filled_df))  

hcc2_dist = dist(agged_mtx, method = "euclidean")
hclust_hcc2 = hclust(hcc2_dist, method = "average")
plot(hclust_hcc2)

agged_mtx <- as.data.frame(agged_mtx)
agged_mtx <- agged_mtx[c("pt3","pt2","pt1"),]
chrs <- paste0('chr', rep(1:22, each=2), '_', rep(c('p', 'q'), 22))
chrs[45] <- "chrX_p"
chrs[46] <- "chrX_q"
# 步骤 4: 生成聚类热图  

pheatmap(combined_df_mtx,
         show_rownames =F,show_colnames = T,
         clustering_method = "mcquitty",cluster_rows = T,cluster_cols = F,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,border_color = 0,
         treeheight_col = 0,angle_col = 315,
         fontsize_col= 12)

combined_df_mtx_t <- t (combined_df_mtx)

row_anno <- as.data.frame(str_split_fixed(row.names(combined_df_mtx_t),"_",2))
row.names(row_anno) <- row.names(combined_df_mtx_t)
row_anno <- as.data.frame(dplyr::select(row_anno,1))
colnames(row_anno) <- c("sample")

col_anno <- as.data.frame(str_split_fixed(colnames(combined_df_mtx_t),"_",2))
row.names(col_anno) <- colnames(combined_df_mtx_t)
col_anno <- as.data.frame(dplyr::select(col_anno,1))
colnames(col_anno) <- c("chr")

ann_colors =list(
  chr=c('chr1'='#E5D2DD','chr2'='#53A85F','chr3'= '#F1BB72','chr4'='#F3B1A0','chr5'='#D6E7A3',
        'chr6'='#57C3F3','chr7'='#476D87','chr8'='#E95C59','chr9'= '#E59CC4','chr10'= '#AB3282', 
        'chr11'='#23452F', 'chr12'='#BD956A', 'chr13'='#8C549C', 'chr14'='#585658', 'chr15'='#9FA3A8',
        'chr16'='#E0D4CA', 'chr17'='#5F3D69','chr18'= '#C5DEBA','chr19'= '#58A4C3','chr20'='#E4C755', 
        'chr21'='#F7F398', 'chr22'='#AA9A59','chrX'= '#E63863', 'chrY'='#E39A35'),
  sample=c('pt1'='#8C549C','pt2'='#D6E7A3','pt3'='#E59CC4')
)

pt3_sub <- subset(row_anno,sample=='pt3')
pt2_sub <- subset(row_anno,sample=='pt2')
pt1_sub <- subset(row_anno,sample=='pt1')

rowname_sort <- c(row.names(pt3_sub),row.names(pt2_sub),row.names(pt1_sub))

pheatmap(combined_df_mtx_t[rowname_sort,],
         show_rownames =F,show_colnames = F,
         scale = 'row',
         clustering_method = "mcquitty",cluster_rows = F,cluster_cols = F,
         annotation_row = row_anno,annotation_col = col_anno,annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,border_color = 0,
         treeheight_col = 0,angle_col = 315,
         fontsize_col= 12)
