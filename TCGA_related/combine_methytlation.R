setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/TCGA_related/methy_merge")
rm(list=ls())

load("methy_merge.RData")

save.image("methy_merge.RData")

trio_meth_mtx2 <- as.data.frame(fread("big_psu_meth.txt",row))
row.names(trio_meth_mtx2) <- trio_meth_mtx2$V1
trio_meth_mtx2 <- trio_meth_mtx2[,-1]
sample_id <- colnames(trio_meth_mtx2)
sample_id <- str_split_i(sample_id,"_D",1)
sample_id <- sample_id[-1]
sample_id <- unique(sample_id)
sample_id <- gsub("hcc3", "HCC2",  sample_id)
sample_id <- gsub("hcc4", "HCC3",  sample_id)
sample_id <- gsub("hcc7", "HCC6",  sample_id)
sample_id <- gsub("hcc11", "HCC7", sample_id)
sample_id <- gsub("hcc28", "HCC8", sample_id)
sample_id <- gsub("hcc29", "HCC9", sample_id)
trio_meth_mtx_sub <- trio_meth_mtx[,sample_id]

trio_meth_mtx <- as.data.frame(fread("trio_meth.txt"))
trio_cites <- row.names(trio_meth_mtx2)
tcga_meth_mtx <- as.data.frame(fread("tcga_meth_100kbin.txt"))
tcga_cites <- tcga_meth_mtx$V1

coverage <- intersect(trio_cites,tcga_cites)
combine_mtx <- cbind(tcga_meth_mtx[coverage,],trio_meth_mtx2[coverage,])
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

trio_col_anno <- as.data.frame(colnames(trio_meth_mtx2))
trio_col_anno[,"origin"] <- "trio-seq"
trio_col_anno[,3:4] <- as.data.frame(str_split_fixed(trio_col_anno[,1],"_",2))
trio_col_anno$V2 <- gsub("\\d+$","Tumor Tissue",trio_col_anno$V2)
trio_col_anno$V2 <- gsub("pt","",trio_col_anno$V2)
trio_col_anno$V2 <- gsub("nt","Normal Tissue",trio_col_anno$V2)
row.names(trio_col_anno) <- trio_col_anno[,1]
trio_col_anno <- trio_col_anno[,-1]
colnames(trio_col_anno) <- c("origin","sample","tissue")



trio_meth_mtx <- as.data.frame(lapply(trio_meth_mtx, as.numeric))

pheatmap(trio_meth_mtx2,
         annotation_col = trio_col_anno,annotation_row = pmd_anno,
         annotation_colors = big_ann_colors_new,
         show_rownames = F,show_colnames = F,
         clustering_method = "mcquitty",
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15)



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
  origin = c("trio-seq"="#defcf9","TCGA" ="#cca8e9")
)
colanno <- select(colanno,!sample)

colanno_geom_point <- colanno

combine_heatmap <-pheatmap(combine_mtx,
         show_rownames = F,show_colnames = F,annotation_col = colanno[,1:2],annotation_row = pmd_anno_psu,
         annotation_colors = big_ann_colors_new,
         clustering_distance_rows = "euclidean",clustering_method = "mcquitty",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 8,
         angle_col = 45)
combine_order_row = combine_heatmap$tree_row$order
combine_order_col = combine_heatmap$tree_col$order
combine_sort_mtx <- data.frame(combine_mtx[combine_order_row, combine_order_col])

tcga_meth_mtx_sub <- tcga_meth_mtx[coverage,]
trio_meth_mtx_sub <- trio_meth_mtx2[coverage,]
sort_trioseq_data = data.frame(trio_meth_mtx_sub[row.names(combine_sort_mtx), colnames(combine_sort_mtx[,431:457])])
sort_tcga_data = data.frame(tcga_meth_mtx_sub[row.names(combine_sort_mtx), colnames(combine_sort_mtx[,1:430])])



pheatmap(sort_trioseq_data,
         show_rownames = F,show_colnames = T,annotation_col = colanno,annotation_row = rowanno,
         annotation_colors = big_ann_colors_new,
         cluster_rows = F,cluster_cols = F,clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 15)











pmd <- fread("PMD_coordinates_hg38.bed")
pmd <- dplyr::select(pmd,1,3,5)
pmd <- mutate(pmd,pos=paste(V1,V3/100000,sep = "_"))
pmd_anno_psu <- as.data.frame(pmd$V5)
colnames(pmd_anno_psu) <- "meth_type"
row.names(pmd_anno_psu) <- pmd$pos
pmd_anno_psu <- na.omit(pmd_anno_psu)

sort_all_id <- colnames(sort_combine_all_data)
sort_tcga_id <- subset(sort_all_id , subset = sort_all_id %in% colnames(tcga_meth_mtx))
sort_tcga <- select(sort_combine_all_data,sort_tcga_id)
sort_trio_id <- subset(sort_all_id , subset = sort_all_id %in% colnames(trio_meth_mtx))
sort_trio <- select(sort_combine_all_data,sort_trio_id)


big_ann_colors_new =list(
  meth_type=c('PMD'='#3f72af','HMD'='#f9ed69'),
  tissue=c('Tumor Tissue'='#b11a2b','Normal Tissue'='#4a74a4'),
  origin=c('TCGA'='#00ADC4','trio-seq'='#4C1A72'),
  mean_methy_level = c("#f7f6ee", "#cca8e9")
)

rownames(colanno)  <- gsub("hcc3", "HCC2",rownames(colanno) )
rownames(colanno) <- gsub("hcc4", "HCC3",rownames(colanno) )
rownames(colanno) <- gsub("hcc7", "HCC6",rownames(colanno) )
rownames(colanno) <- gsub("hcc11", "HCC7",rownames(colanno) )
rownames(colanno) <- gsub("hcc28", "HCC8",rownames(colanno) )
rownames(colanno) <- gsub("hcc29", "HCC9",rownames(colanno)  )
colnames(colanno) <- c("origin","tissue","mean")
colanno$tissue <- gsub("tumor", "Tumor Tissue",colanno$tissue)
colanno$tissue <- gsub("normal", "Normal Tissue",colanno$tissue)
colanno['HCC6_pt1',] <- c("trio-seq","Tumor Tissue",0.7)


pmd_regions <- subset(pmd_anno_psu,subset=meth_type =='PMD') 
pmd_regions$region <- row.names(pmd_regions)
pmd_regions <-subset(pmd_regions ,subset=region %in% row.names(sort_tcga_data) )

hmd_regions <- subset(pmd_anno_psu,subset=meth_type =='HMD') 
hmd_regions$region <- row.names(hmd_regions)
hmd_regions <-subset(hmd_regions ,subset=region %in% row.names(sort_tcga_data) )


sorted_tcga_colmeans <- as.data.frame(colMeans(sort_tcga_data))
colnames(sorted_tcga_colmeans) <- 'mean_level'
sorted_tcga_colmeans$pmd_level <- sorted_tcga_pmd_colmeans$`colMeans(sort_tcga_data[rownames(pmd_regions), ])`
sorted_tcga_colmeans$hmd_level <- sorted_tcga_hmd_colmeans$hmd_level

sorted_trioseq_colmeans <- as.data.frame(colMeans(sort_trioseq_data))
colnames(sorted_trioseq_colmeans) <- 'mean_level'
sorted_trioseq_colmeans$pmd_level <- sorted_trioseq_pmd_colmeans$`colMeans(sort_trioseq_data[rownames(pmd_regions), ])`
sorted_trioseq_colmeans$hmd_level <- sorted_trioseq_hmd_colmeans$hmd_level

sorted_tcga_pmd_colmeans <-as.data.frame(colMeans(sort_tcga_data[rownames(pmd_regions),]))
colnames(sorted_tcga_pmd_colmeans) <- 'pmd_level'
sorted_trioseq_pmd_colmeans <-as.data.frame(colMeans(sort_trioseq_data[rownames(pmd_regions),]))
colnames(sorted_trioseq_pmd_colmeans) <- 'pmd_level'

sorted_tcga_hmd_colmeans <-as.data.frame(colMeans(sort_tcga_data[rownames(hmd_regions),]))
colnames(sorted_tcga_hmd_colmeans) <- 'hmd_level'
sorted_trioseq_hmd_colmeans <-as.data.frame(colMeans(sort_trioseq_data[rownames(hmd_regions),]))
colnames(sorted_trioseq_hmd_colmeans) <- 'hmd_level'

meanlevel_sum <- rbind(sorted_trioseq_colmeans,sorted_tcga_colmeans)
meanlevel_sum$sample <- row.names(meanlevel_sum) 

meanlevel_sum$sample <-factor(meanlevel_sum$sample, levels =c(meanlevel_sum$sample))




ggplot(meanlevel_sum[34:457,])+
  geom_point(aes(x=sample,y=hmd_level),color='#f9ed69',size=1)+
  geom_line(aes(x=factor(sample),y=hmd_level,group=1),size = 1,color='#f9ed69')+
  geom_line(aes(x=factor(sample),y=pmd_level,group=1),size = 1,color='#3f72af')+
  geom_point(aes(x=sample,y=pmd_level),color='#3f72af',size=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.title.x = element_text(size = 0),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 0,color="black",angle = 45),axis.text.y = element_text(size = 10,color="black"))


ggplot(meanlevel_sum)+
  geom_point(aes(x=sample,y=hmd_level),color='#f9ed69',size=1)+
  geom_line(aes(x=factor(sample),y=hmd_level,group=1),size = 1,color='#f9ed69')+
  geom_line(aes(x=factor(sample),y=pmd_level,group=1),size = 1,color='#3f72af')+
  geom_point(aes(x=sample,y=pmd_level),color='#3f72af',size=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.title.x = element_text(size = 0),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 0,color="black",angle = 45),axis.text.y = element_text(size = 10,color="black"))




colanno$mean <- as.numeric(colanno$mean)

colnames(colanno) <- c("origin","tissue","mean_methy_level")
cg_count <- fread("cg_count.txt")
cg_count <- subset(cg_count,subset=Var1 %in% row.names(sort_tcga_data))

cg_count$Var1 <- factor(cg_count$Var1,levels = row.names(sort_tcga_data))

ggplot(cg_count)+
  geom_line(aes(x=factor(Var1),y=Freq,group=1),size = 0.4,color='black')

pheatmap(sort_tcga_data,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = pmd_anno_psu,cluster_cols = F,cluster_rows = F,
         clustering_method = "mcquitty",annotation_colors = big_ann_colors_new,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 1,
         cellheight = 0.03)
pheatmap(sort_trioseq_data,
         show_rownames = F,show_colnames = F,annotation_col = colanno,annotation_row = pmd_anno_psu,cluster_cols = F,cluster_rows = F,
         clustering_method = "mcquitty",annotation_colors = big_ann_colors_new,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 10,
         cellheight = 0.03)

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
all_demeth <- demeth1


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




all_demeth_regions <- as.data.frame(str_split_fixed(row.names(all_demeth),"_",2))
colnames(all_demeth_regions) <- c("chr","start")
all_demeth_regions$pos <- row.names(all_demeth)
all_demeth_regions$start <- as.numeric(all_demeth_regions$start)
all_demeth_regions$end <- all_demeth_regions$start+1
all_demeth_regions$start <- all_demeth_regions$start*100000
all_demeth_regions$end <- all_demeth_regions$end*100000
all_demeth_regions$end <- as.numeric(all_demeth_regions$end)
all_demeth_regions$start <- as.numeric(all_demeth_regions$start)
all_demeth_regions <- select(all_demeth_regions,chr,start,end,pos)
write.table(all_demeth_regions, "all_demeth_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)





trio_meth_10x <- select(trio_meth_mtx,hcc3_nt,hcc3_pt1,hcc3_pt2,hcc3_pt3,hcc3_pt4,hcc28_nt,hcc28_pt1,hcc28_pt2,hcc28_pt4,hcc29_pt1,hcc29_pt3,hcc29_pt4)
pheatmap(trio_meth_10x,
         show_rownames = F,show_colnames = T,annotation_row = rowanno,
         annotation_colors = merge_ann_colors_sub,cluster_cols = F,
         clustering_distance_rows = "euclidean",clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 40,
         fontsize_col= 13,
         angle_col = 45,
         cellwidth = 25)


trio_meth_hcc3 <- select(sort_trioseq_data,HCC2_pt1,HCC2_pt2,HCC2_pt3)
trio_meth_hcc3_pt <- select(trio_meth_mtx,hcc3_pt1,hcc3_pt2,hcc3_pt3,hcc3_pt4)
pheatmap(trio_meth_hcc3,
         show_rownames = F,show_colnames = T,annotation_row = pmd_anno_psu,
         annotation_colors = big_ann_colors_new,cluster_cols = F,
         clustering_distance_rows = "euclidean",clustering_method = "average",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 40,
         fontsize_col= 13,
         angle_col = 45,
         cellwidth = 80)

pheatmap(trio_meth_hcc3_pt,
         show_rownames = F,show_colnames = T,
         annotation_row = pmd_anno_psu,
         annotation_colors = merge_ann_colors_sub,
         clustering_distance_rows = "euclidean",clustering_method = "average",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 100,
         fontsize_col= 13,
         angle_col = 45,
         cellwidth = 80)
hcc3heatmap_result <- pheatmap(trio_meth_hcc3_pt,
                               show_rownames = F,show_colnames = T,annotation_row = rowanno,
                               annotation_colors = merge_ann_colors_sub,
                               clustering_distance_rows = "euclidean",clustering_method = "average",
                               color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                               treeheight_row = 0,
                               treeheight_col = 100,
                               fontsize_col= 13,
                               angle_col = 45,
                               cellwidth = 80)
hcc3_pb_mtx <- t(trio_meth_hcc3)
hcc3_pb_dist = dist(hcc3_pb_mtx, method = "euclidean")
hclust_hcc3 = hclust(hcc3_pb_dist, method = "average")
plot(hclust_hcc3)

plot(hhcc3_distance_matrix)



trio_meth_hcc28 <- select(trio_meth_mtx,hcc28_nt,hcc28_pt1,hcc28_pt2,hcc28_pt4)

pheatmap(trio_meth_hcc28,
         show_rownames = F,show_colnames = T,annotation_row = rowanno,
         annotation_colors = merge_ann_colors_sub,
         clustering_distance_rows = "euclidean",clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 15,
         fontsize_col= 13,
         angle_col = 45,
         cellwidth = 80)



trio_meth_hcc28 <- select(trio_meth_mtx,hcc28_pt1,hcc28_pt2,hcc28_pt4)

pheatmap(trio_meth_hcc28,
         show_rownames = F,show_colnames = T,annotation_row = rowanno,
         annotation_colors = merge_ann_colors_sub,
         clustering_distance_rows = "euclidean",clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 15,
         fontsize_col= 13,
         angle_col = 45,
         cellwidth = 80)

hcc28_pb_mtx <- t(trio_meth_hcc28)
hcc28_pb_dist = dist(hcc28_pb_mtx, method = "euclidean")
hclust_hcc28 = hclust(hcc28_pb_dist, method = "average")
plot(hclust_hcc28)






trio_meth_hcc29 <- select(trio_meth_mtx,hcc29_pt1,hcc29_pt3,hcc29_pt4)

pheatmap(trio_meth_hcc29,
         show_rownames = F,show_colnames = T,annotation_row = rowanno,
         annotation_colors = merge_ann_colors_sub,
         clustering_distance_rows = "euclidean",clustering_method = "average",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 15,
         fontsize_col= 13,
         angle_col = 45,
         cellwidth = 80)

hcc29_pb_mtx <- t(trio_meth_hcc29)
hcc29_pb_dist = dist(hcc29_pb_mtx, method = "euclidean")
hclust_hcc29 = hclust(hcc29_pb_dist, method = "average")
plot(hclust_hcc29)


setwd("~/projects/hcc/data/trio_seq/scMethy/methy_data/methy_file/hcc7/CpG_profile/pt1")
hcc6_pt1_d1 <- fread("pt1_D1.CpG.methy.bed.gz")
setwd("~/projects/hcc/data/trio_seq/scMethy/methy_data/methy_file/hcc7/CpG_100kb")
hcc6_pt1_d1 <- fread("pt1_D1.CpG.100kb.methy.tsv")
