cor_matrix

bigmeth_cor <- big_meth

colnames(bigmeth_cor) <- gsub("hcc3", "HCC2",colnames(bigmeth_cor) )
colnames(bigmeth_cor) <- gsub("hcc4", "HCC3",colnames(bigmeth_cor) )
colnames(bigmeth_cor) <- gsub("hcc7", "HCC6",colnames(bigmeth_cor) )
colnames(bigmeth_cor) <- gsub("hcc11", "HCC7",colnames(bigmeth_cor) )
colnames(bigmeth_cor) <- gsub("hcc28", "HCC8",colnames(bigmeth_cor) )
colnames(bigmeth_cor) <- gsub("hcc29", "HCC9",colnames(bigmeth_cor)  )
cor_list <- list()
sample_names_cor <- sub("^(([^_]+_[^_]+)).*", "\\1", colnames(bigmeth_cor))
unique_samples <- unique(sample_names_cor)
patient_id <- sub("_.*", "", colnames(bigmeth_cor))


for (s in unique_samples) {
  cat("Processing sample:", s, "\n")
  
  # 提取该样本对应的列
  cells_in_sample <- colnames(bigmeth_cor)[sample_names_cor == s]
  if (length(cells_in_sample) < 2) {
    cat("Skip sample", s, "because it has less than 2 cells\n")
    next
  }
  
  sub_matrix <- bigmeth_cor[, cells_in_sample, drop=FALSE]
  
  cor_mat <- cor(sub_matrix, use="pairwise.complete.obs", method="pearson")
  cor_values <- cor_mat[upper.tri(cor_mat)]
  
  # 提取病人编号
  pid <- sub("_.*", "", s)
  
  cor_list[[s]] <- data.frame(
    sample = s,
    patient = pid,
    cor_value = cor_values
  )
}

cor_df <- bind_rows(cor_list)

cor_df <- mutate(cor_df)

patient_type <- data.frame(
  patient = c("HCC2", "HCC6", "HCC8", "HCC9", "HCC3", "HCC7"),
  type = c("CMN", "CMN", "CMN", "CMN", "SN", "SN")
)


cor_df <- left_join(cor_df, patient_type, by="patient")
library(dplyr)

cor_df <- cor_df %>%
  mutate(group = ifelse(sub("^[^_]+_([a-zA-Z]).*", "\\1", sample) == "n", "normal", "tumor"))

sta_samples <- c("HCC3_pt3", "HCC6_pt4", "HCC6_pt5", "HCC6_pt6", "HCC7_pt6")

cor_df <- cor_df %>%
  mutate(pos = ifelse(sample %in% sta_samples, "sta", "prm"))



cor_df_t <- subset(cor_df,subset=group=="tumor")
ggplot(cor_df_t,aes(type,cor_value,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("SN","CMN")),
                     method = "wilcox.test",label = "p.format",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

cor_df_st <- subset(cor_df_t,subset=patient %in% c("HCC3","HCC6","HCC7"))

ggplot(cor_df_st,aes(pos,cor_value,color=pos))+
  geom_boxplot(width=0.5,outlier.size=0)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("sta","prm")),
                     method = "wilcox.test",label = "p.format",
                     label.y =1 )+theme_bw() +facet_wrap(~ patient, scales = "free_y") 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())
  
  
  
ggplot(cor_df_st,aes(type,cor_value,color=type))+
    geom_boxplot(width=0.5,outlier.size=0)+
    scale_color_manual(values =c('#b11a2b','#4a74a4'))+
    stat_compare_means(comparisons = list(c("CMN","SN")),
                       method = "wilcox.test",label = "p.format",
                       label.y =1 )+theme_bw() +facet_wrap(~ pos, scales = "free_y")
  
  
