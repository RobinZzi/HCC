setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/TCGA_related/methy_merge")
rm(list=ls())

load("methy_merge.RData")

save.image("methy_merge.RData")


trio_meth_mtx <- as.data.frame(fread("trio_meth.txt"))
trio_cites <- trio_meth_mtx$V1
tcga_meth_mtx <- as.data.frame(fread("tcga_meth_100kbin.txt"))
tcga_cites <- tcga_meth_mtx$V1

coverage <- intersect(trio_meth_mtx$V1,tcga_meth_mtx$V1)
combine_mtx <- merge(tcga_meth_mtx,trio_meth_mtx,by='V1')
coverage_cites <- combine_mtx$V1

trio_meth_mtx <- trio_meth_mtx[,-1]
trio_meth_mtx <- as.data.frame(lapply(trio_meth_mtx,as.numeric))
row.names(trio_meth_mtx) <- trio_cites

tcga_meth_mtx <- tcga_meth_mtx[,-1]
tcga_meth_mtx <- as.data.frame(lapply(tcga_meth_mtx,as.numeric))
row.names(tcga_meth_mtx) <- tcga_cites

combine_mtx <- combine_mtx[,-1]
combine_mtx <- as.data.frame(lapply(combine_mtx,as.numeric))
row.names(combine_mtx) <- coverage_cites


tcga_col_anno <- as.data.frame(colnames(tcga_meth_mtx))
tcga_col_anno[,"origin"] <- "TCGA"
tcga_col_anno[,"sample_type"] <- str_sub(tcga_col_anno$`colnames(tcga_meth_mtx)`,14,15)
tcga_col_anno[,"sample_type"] <- gsub("01","tumor",tcga_col_anno[,"sample_type"])
tcga_col_anno[,"sample_type"] <- gsub("02","tumor",tcga_col_anno[,"sample_type"])
tcga_col_anno[,"sample_type"] <- gsub("11","normal",tcga_col_anno[,"sample_type"])
tcga_col_anno <- tcga_col_anno[-1,]
row.names(tcga_col_anno) <- tcga_col_anno[,1]
tcga_col_anno <- tcga_col_anno[,-1]
tcga_col_anno$sample <- "TCGA"
tcga_col_anno <- select(tcga_col_anno,origin,sample,sample_type)
row.names(tcga_col_anno) <- gsub("-",'.',row.names(tcga_col_anno)) 

trio_col_anno <- as.data.frame(colnames(trio_meth_mtx))
trio_col_anno[,"origin"] <- "trio-seq"
trio_col_anno[,3:4] <- as.data.frame(str_split_fixed(trio_col_anno$`colnames(trio_meth_mtx)`,"_",2))
trio_col_anno$V2 <- gsub("\\d+$","tumor",trio_col_anno$V2)
trio_col_anno$V2 <- gsub("pt","",trio_col_anno$V2)
trio_col_anno$V2 <- gsub("nt","normal",trio_col_anno$V2)
row.names(trio_col_anno) <- trio_col_anno[,1]
trio_col_anno <- trio_col_anno[,-1]
colnames(trio_col_anno) <- c("origin","sample","sample_type")














tcga_cover_graph <- as.data.frame(str_split_fixed(row.names(tcga_meth_mtx),"_",2))
colnames(tcga_cover_graph) <- c("chr","start")
tcga_cover_graph$pos <- row.names(tcga_meth_mtx)
tcga_cover_graph$start <- as.numeric(tcga_cover_graph$start)
tcga_cover_graph$end <- tcga_cover_graph$start+1
tcga_cover_graph$start <- tcga_cover_graph$start*100000
tcga_cover_graph$end <- tcga_cover_graph$end*100000
tcga_cover_graph$end <- as.numeric(tcga_cover_graph$end)
tcga_cover_graph$start <- as.numeric(tcga_cover_graph$start)
tcga_cover_graph <- select(tcga_cover_graph,chr,start,end,pos)
write.table(tcga_cover_graph, "tcga_cover.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)







pwd_anno <- fread("tcga_PMD.bedGraph")
pwd_anno  <- select(pwd_anno ,V5,V6,V10)
pwd_anno  <- aggregate(pwd_anno,by=list(pwd_anno$V10),max)
row.names(pwd_anno) <- pwd_anno$Group.1 
pwd_anno <- select(pwd_anno,V5,V6)
colnames(pwd_anno) <- c("meth_type","pmd_type")
pmd_anno <- pwd_anno
pmd_anno <- select(pmd_anno,meth_type)

tcga_pmd_region <- subset(pmd_anno,subset = meth_type=="PMD")
pmd_sub_mtx <- tcga_meth_mtx[row.names(tcga_pmd_region),]
sub_colanno <- as.data.frame(colMeans(pmd_sub_mtx))
colnames(sub_colanno) <- "PMD_meanmeth_level"

sub_colanno <- merge(tcga_col_anno,sub_colanno,by="row.names")
row.names(sub_colanno) <- sub_colanno$Row.names
sub_colanno <-  sub_colanno[,-1]
sub_colanno <- select(sub_colanno,sample_type,PMD_meanmeth_level)
pheatmap(tcga_meth_mtx,
         show_rownames = F,show_colnames = F,annotation_col = sub_colanno,annotation_row = pmd_anno,
         annotation_colors = merge_ann_colors_sub,
         cluster_rows = T,cluster_cols = T,clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 1)



colanno <- rbind(tcga_col_anno,trio_col_anno)
combine_pmd_region <- subset(rowanno,subset = meth_type=="PMD")
pmd_combine_mtx <- combine_mtx[row.names(combine_pmd_region),]
pmd_combine_mtx <- na.omit(pmd_combine_mtx)
combine_meanlevel <- as.data.frame(colMeans(pmd_combine_mtx))
colnames(combine_meanlevel) <- "PMD_meanmeth_level"
colanno <- merge(colanno,combine_meanlevel,by="row.names")
rownames(colanno) <- colanno$Row.names
colanno <- colanno[,-1]

rowanno <- as.data.frame(fread("row_anno.txt"))
row.names(rowanno) <- rowanno$V1
rowanno <- rowanno[,-1]
rowanno <- select(rowanno,meth_type)
merge_ann_colors_sub=list(
  rep_time_norm = colorRampPalette(c("#28DCDC", "#fcefee", "#d72323"))(200),
  rep_time=c('LRD'='#99CCFF','ERD'='#FF0000','DTZ'='#FFE0EC','UTZ'='#FFF599'),
  pmd_type=c('commonPMD'='#99CCFF','commonHMD'='#FF0000','Neither'='#FFE0EC'),
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  sample_type = c("normal"='#99CCFF',"tumor"='#FF0000'),
  origin = c("trio-seq"="#C8D948","TCGA" ="#84C7DB")
)
colanno <- select(colanno,!sample)

colanno_geom_point <- colanno

pheatmap(combine_mtx,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
         annotation_colors = merge_ann_colors_sub,
         clustering_distance_rows = "euclidean",clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 1)
sub_trio_meth_mtx <- trio_meth_mtx[coverage_cites,]
trioseq_heatmap <- pheatmap(sub_trio_meth_mtx,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
         clustering_method = "mcquitty",annotation_colors = merge_ann_colors_sub,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 15)

trioseq_order_row = trioseq_heatmap$tree_row$order
trioseq_order_col = trioseq_heatmap$tree_col$order
sort_trioseq_data = data.frame(sub_trio_meth_mtx[trioseq_order_row, trioseq_order_col])
sort_combine_data = data.frame(combine_mtx[trioseq_order_row, ])


pheatmap(sort_trioseq_data,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
         annotation_colors = merge_ann_colors_sub,
         cluster_rows = F,cluster_cols = F,clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 15)



pheatmap(sort_combine_data,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
         annotation_colors = merge_ann_colors_sub,
         cluster_rows = F,cluster_cols = T,clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 1)


trioseq_all_heatmap <- pheatmap(trio_meth_mtx,
                            show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
                            clustering_method = "mcquitty",annotation_colors = merge_ann_colors_sub,
                            clustering_distance_rows = "euclidean",
                            color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                            treeheight_row = 0,
                            treeheight_col = 1,
                            fontsize_col= 8,
                            angle_col = 45,
                            cellwidth = 15)
trioseq_all_order_row = trioseq_all_heatmap$tree_row$order
trioseq_all_order_col = trioseq_all_heatmap$tree_col$order
sort_trioseq_all_data = data.frame(trio_meth_mtx[trioseq_all_order_row, trioseq_all_order_col])

dev.off()

combine_all_heatmap <- pheatmap(combine_mtx,
                                show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
                                clustering_method = "mcquitty",annotation_colors = merge_ann_colors_sub,
                                clustering_distance_rows = "euclidean",
                                color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                                treeheight_row = 0,
                                treeheight_col = 1,
                                fontsize_col= 8,
                                angle_col = 45,
                                cellwidth = 1)
combine_all_order_row = combine_all_heatmap$tree_row$order
combine_all_order_col = combine_all_heatmap$tree_col$order
sort_combine_all_data = data.frame(combine_mtx[combine_all_order_row, combine_all_order_col])

colanno_geom_point <- colanno
colanno_geom_point$sample <- row.names(colanno_geom_point)

ggplot(colanno_geom_point,aes(x=factor(sample,levels = colnames(sort_combine_all_data)),y=PMD_meanmeth_level))+
  geom_point(size = 0.5)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "", 
       y= " ",
       title="")


tcga_all_heatmap <- pheatmap(tcga_meth_mtx,
                                show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
                                clustering_method = "mcquitty",annotation_colors = merge_ann_colors_sub,
                                clustering_distance_rows = "euclidean",
                                color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                                treeheight_row = 0,
                                treeheight_col = 1,
                                fontsize_col= 8,
                                angle_col = 45,
                                cellwidth = 1)
tcga_order_row = tcga_all_heatmap$tree_row$order
tcga_order_col = tcga_all_heatmap$tree_col$order
sort_tcga_all_data = data.frame(tcga_meth_mtx[tcga_order_row, tcga_order_col])

colanno_geom_point_tcga <- subset (colanno_geom_point,subset = origin == 'TCGA')


ggplot(colanno_geom_point_tcga,aes(x=factor(sample,levels = colnames(sort_tcga_all_data)),y=PMD_meanmeth_level))+
  geom_point(size = 0.5)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "", 
       y= " ",
       title="")

demeth1 <- sort_trioseq_all_data[6835:16747,]
pheatmap(demeth1,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
         cluster_rows = T,cluster_cols = T,annotation_colors = merge_ann_colors_sub,
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 15)















demeth1 <- sort_trioseq_all_data[6835:16747,]
pheatmap(demeth1,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = rowanno,
         cluster_rows = T,cluster_cols = T,annotation_colors = merge_ann_colors_sub,
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 15)


late_demeth <- demeth1[1:5824,]
early_demeth <- demeth1[5825:9913,]



late_demeth_regions <- as.data.frame(str_split_fixed(row.names(late_demeth),"_",2))
colnames(late_demeth_regions) <- c("chr","start")
late_demeth_regions$pos <- row.names(late_demeth)
late_demeth_regions$start <- as.numeric(late_demeth_regions$start)
late_demeth_regions$end <- late_demeth_regions$start+1
late_demeth_regions$start <- late_demeth_regions$start*100000
late_demeth_regions$end <- late_demeth_regions$end*100000
late_demeth_regions$end <- as.numeric(late_demeth_regions$end)
late_demeth_regions$start <- as.numeric(late_demeth_regions$start)
late_demeth_regions <- select(late_demeth_regions,chr,start,end,pos)
write.table(late_demeth_regions, "late_demeth_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



early_demeth_regions <- as.data.frame(str_split_fixed(row.names(early_demeth),"_",2))
colnames(early_demeth_regions) <- c("chr","start")
early_demeth_regions$pos <- row.names(early_demeth)
early_demeth_regions$start <- as.numeric(early_demeth_regions$start)
early_demeth_regions$end <- early_demeth_regions$start+1
early_demeth_regions$start <- early_demeth_regions$start*100000
early_demeth_regions$end <- early_demeth_regions$end*100000
early_demeth_regions$end <- as.numeric(early_demeth_regions$end)
early_demeth_regions$start <- as.numeric(early_demeth_regions$start)
early_demeth_regions <- select(early_demeth_regions,chr,start,end,pos)
write.table(early_demeth_regions, "early_demeth_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
