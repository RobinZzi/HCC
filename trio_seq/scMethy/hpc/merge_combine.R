setwd("/storage/zhangyanxiaoLab/zhangliwen/Projects/sc-LIHC/analysis/scMehty/")
library(data.table)
library(stringr)
library(dplyr)
library(pheatmap)
rm(list=ls())
fs <- list.files('./',pattern = "*.txt")
save.image("merge.Rdata")
for(i in 1:length(fs)){
  data <- fread(fs[i])
  id <- str_split(paste(fs[i]),".txt")[[1]][1]
  colnames(data) <- c("chr","site","meth","dmeth")
  data[,1:2] <-data.frame( str_split_fixed(data$site,"_",2))
  data$site <- as.numeric(data$site)
  data[,2] <- round(data[,2]/100000)
  data <- tidyr::unite(data,"id",chr,site,sep="_",remove=TRUE)
  data <- aggregate(data[,2:3],by=list(id=data$id),FUN=sum)
  data[,paste(id)] <- data[,2]/(data[,2]+data[,3])
  data <- select(data,id,paste(id))
  if(i == 1){
    base <- data
  }else{
    base <- merge(base,data,by="id")
  }
  print(paste(i,id))
}
write.table(base,"merge.txt")






merge_result <- fread("merge.txt")
merge_result <- merge_result[,-1]
merge_result <- as.data.frame(merge_result)
row.names(merge_result) <- merge_result$id
merge_result <- merge_result[,-1]
merge_result <- select(merge_result,-hcc7_pt6)
























merge_graph <- as.data.frame(str_split_fixed(row.names(merge_result),"_",2))
colnames(merge_graph) <- c("chr","start")
merge_graph$pos <- row.names(merge_result)
merge_graph$start <- as.numeric(merge_graph$start)
merge_graph$end <- merge_graph$start+1
merge_graph$start <- merge_graph$start*100000
merge_graph$end <- merge_graph$end*100000
merge_graph$end <- as.numeric(merge_graph$end)
merge_graph$start <- as.numeric(merge_graph$start)
merge_graph <- select(merge_graph,chr,start,end,pos)
write.table(merge_graph, "merge.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)









merge_rowanno2 <- fread("merge_time2.bedGraph")
merge_rowanno2 <- select(merge_rowanno2,V4,V8,V9)
merge_rowanno2 <- aggregate(merge_rowanno2,by=list(merge_rowanno2$V8),max)
row.names(merge_rowanno2) <- merge_rowanno2$Group.1 
merge_rowanno2 <- select(merge_rowanno2,V4)
colnames(merge_rowanno2) <- "rep_time"
merge_rowanno2 <- as.data.frame(merge_rowanno2)


merge_rowanno2$rep_time_norm <- scale(merge_rowanno2$rep_time)
merge_rowanno2$rep_time_log <- log10(merge_rowanno2$rep_time)


merge_rowanno2$rep_time_norm


merge_rowanno3$rep_time_norm[merge_rowanno3$rep_time_norm >= 2] = 2
merge_rowanno3$rep_time_norm[merge_rowanno3$rep_time_norm <= -2] = -2





colnames(merge_rowanno2) <- c("rep_time","rep_time_norm","rep_time_log")

merge_rowanno <- fread("merge_time.bedGraph")
merge_rowanno <- select(merge_rowanno,V4,V8,V9)
merge_rowanno <- aggregate(merge_rowanno,by=list(merge_rowanno$V8),max)
row.names(merge_rowanno) <- merge_rowanno$Group.1 
merge_rowanno <- select(merge_rowanno,V4)
colnames(merge_rowanno) <- "rep_time"


merge_ann_colors=list(
  rep_time=c('LRD'='#0039B3','ERD'='#FF1962','DTZ'='#FFE0EC','UTZ'='#FFF599')
)

merge_rowanno$rep_time <-scale(merge_rowanno$rep_time, center=T,scale=T)
colnames(merge_rowanno) <- "rep_time"

merge_ann_colors2=list(
  rep_time = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(200),
  rep_time_log = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(200),
  rep_time_norm = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(200)
)

merge_rowanno3 <- select(merge_rowanno2,rep_time_norm)
merge_ann_colors3=list(
  rep_time_norm = colorRampPalette(c("#2EB3B0", "#fcefee", "#d72323"))(200)
)

pheatmap(merge_result,
         show_rownames = F,show_colnames = T,
         clustering_method = "mcquitty",annotation_row = merge_rowanno3,
         clustering_distance_rows = "euclidean",annotation_colors = merge_ann_colors3,
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 12)











solo_pos_anno <- fread("merge_solo.bedGraph")
solo_pos_anno <- select(solo_pos_anno,V4,V8,V9)
solo_pos_anno <- aggregate(solo_pos_anno,by=list(solo_pos_anno$V9),max)
row.names(solo_pos_anno) <- solo_pos_anno$Group.1 
solo_pos_anno <- select(solo_pos_anno,V4)
colnames(solo_pos_anno) <- "solo_pos"


pwd_anno <- fread("merge_pwd.bedGraph")
pwd_anno  <- select(pwd_anno ,V5,V6,V10)
pwd_anno  <- aggregate(pwd_anno,by=list(pwd_anno$V10),max)
row.names(pwd_anno) <- pwd_anno$Group.1 
pwd_anno <- select(pwd_anno,V5,V6)
colnames(pwd_anno) <- c("meth_type","pmd_type")






pheatmap(merge_result,
         show_rownames = F,show_colnames = T,
         clustering_method = "mcquitty",annotation_row = row_anno_sub,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 12)




row_anno_sub <- merge(pwd_anno,merge_rowanno3,by= "row.names")
rownames(row_anno_sub) <- row_anno_sub$Row.names
row_anno_sub <- row_anno_sub[,-1]
row_anno_sub <- merge(row_anno_sub,merge_rowanno,by= "row.names")
rownames(row_anno_sub) <- row_anno_sub$Row.names
row_anno_sub <- row_anno_sub[,-1]




merge_ann_colors_sub=list(
  rep_time_norm = colorRampPalette(c("#28DCDC", "#fcefee", "#d72323"))(200),
  rep_time=c('LRD'='#99CCFF','ERD'='#FF0000','DTZ'='#FFE0EC','UTZ'='#FFF599'),
  pmd_type=c('commonPMD'='#99CCFF','commonHMD'='#FF0000','Neither'='#FFE0EC'),
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000')
)
colanno <- as.data.frame(colMeans(merge_result))
colnames(colanno) <- "mean_level"



pheatmap(merge_result,
         show_rownames = F,show_colnames = T,
         clustering_method = "mcquitty",annotation_row = row_anno_sub,annotation_col = colanno,
         clustering_distance_rows = "euclidean",annotation_colors = merge_ann_colors_sub,
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 15)
heatmap <- pheatmap(merge_result,
                    show_rownames = F,show_colnames = T,
                    clustering_method = "mcquitty",annotation_row = row_anno_sub,annotation_col = colanno,
                    clustering_distance_rows = "euclidean",annotation_colors = merge_ann_colors_sub,
                    color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                    treeheight_row = 0,
                    treeheight_col = 1,
                    fontsize_col= 8,
                    angle_col = 45,
                    cellwidth = 15)
