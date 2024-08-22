
setwd("~/projects/hcc/analysis/inferCNV/clustering")

rm(list = ls())
save.image("hcc2_infer.Rdata")

#data-load#
setwd("~/projects/hcc/analysis/inferCNV/hpc_rds/10x")
hcc2_hpc <- readRDS("~/projects/hcc/analysis/inferCNV/hpc_rds/10x/hcc2_hpc.RDS")
setwd("~/projects/hcc/analysis/inferCNV/infer_result/10x/hcc2_hpc")
hcc2_mtx <- as.data.frame(fread("final.infercnv_matrix.csv"))
row.names(hcc2_mtx) <- hcc2_mtx$V1
hcc2_mtx <- hcc2_mtx[,-1]

hcc2_meta <- as.data.frame(hcc2_hpc$sample_pt)
colnames(hcc2_meta) <- 'sample_pt'
hcc2_meta$Cell_ID <- row.names(hcc2_meta)





# 转换基因表达矩阵为长格式  
long_expression_data <- hcc2_mtx %>%  
  as.data.frame() %>%  
  rownames_to_column("Gene") %>%  
  pivot_longer(cols = -Gene, names_to = "Cell_ID", values_to = "Expression")  

# 合并数据框，加入样本信息  
merged_expression_data <- long_expression_data %>%  
  left_join(hcc2_meta, by = "Cell_ID")  

# 计算每个样本每个基因的平均表达水平  
average_expression_by_sample <- merged_expression_data %>%  
  group_by(sample_pt, Gene) %>%  
  summarize(Average_Expression = mean(Expression, na.rm = TRUE), .groups = 'drop')  

# 转换回宽格式  
wide_average_expression <- average_expression_by_sample %>%  
  pivot_wider(names_from = sample_pt, values_from = Average_Expression)  

# 将基因设为行名  
wide_average_expression <- column_to_rownames(wide_average_expression, "Gene")

hcc2_average_cnv_by_sample_sub <- hcc2_average_cnv_by_sample[c('hcc3_PT1','hcc3_PT2','hcc3_PT3'),]

# 转置以适应距离矩阵计算  
hcc2_average_cnv_by_sample <- t(wide_average_expression)  
hcc2_average_cnv_by_sample <- as.matrix(hcc2_average_cnv_by_sample)  
# 计算样本之间的距离矩阵  
hcc2_distance_matrix <- dist(hcc2_average_cnv_by_sample, method="maximum")
kmeans_hcc2 <- kmeans(hcc2_average_cnv_by_sample_sub)
hclust_hcc2 = hclust(hcc2_distance_matrix, method = "average")
plot(hclust_hcc2)
pheatmap(hcc2_average_cnv_by_sample_sub,show_rownames = T,show_colnames = F)
