setwd("~/projects/hcc/data/wes_hcc/hcc3/gatk_result/cnvkit_result")

hcc3_PT1 <- fread("BCPT1_bqsr.cns")
hcc3_PT2 <- fread("BCPT2_bqsr.cns")
hcc3_PT3 <- fread("BCPT3_bqsr.cns")
hcc3_PT4 <- fread("BCPT4_bqsr.cns")



library(dplyr)  
library(readr)  
library(tidyverse)
library(pheatmap)  

# 步骤 1: 读取所有 cnr 文件  
path_to_cnr_files <- getwd() 
file_list <- list.files(path_to_cnr_files, pattern="*.cnr", full.names=TRUE)  

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

# 步骤 2: 移除不需要的列，仅保留 log2 ratios  
combined_df <- combined_df %>% select(-chromosome, -start, -end)  

# 步骤 3: 填补缺失值并标准化数据  
filled_df <- combined_df %>% replace(is.na(.), 0)  
standardized_df <- scale(t(filled_df))  
standardized_df_t <- t(standardized_df)
hcc3_wes_dist = dist(standardized_df_t, method = "euclidean")
hclust_hcc3 = hclust(hcc3_wes_dist, method = "average")

# 步骤 4: 生成聚类热图  
pheatmap(standardized_df,   
         clustering_method = "ward.D2",   
         clustering_distance_rows = "euclidean",   
         clustering_distance_cols = "euclidean",   
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
