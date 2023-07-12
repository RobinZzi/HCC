rm(list=ls())
library(xlsx)
library(stringr)
library(pheatmap)
library(Seurat)
library(ggforce)
library(reshape2)

library(geneRal)
setwd("~/projects/hcc/analysis/trio_seq/scMethy")
rm(list=ls())
load("methy_merge.Rdata")
save.image("methy_merge.Rdata")



hcc_3_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc3") 
hcc_4_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc4")
hcc_7_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc7") 
hcc_11_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc11") 
hcc_28_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc28") 
hcc_29_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc29") 


hcc_3_qc_pass <- subset(hcc_3_qc, subset= hcc_3_qc$最终是否可用 ==1)
hcc_4_qc_pass <- subset(hcc_4_qc, subset= hcc_4_qc$最终是否可用 ==1)
hcc_7_qc_pass <- subset(hcc_7_qc, subset= hcc_7_qc$最终是否可用 ==1)
hcc_11_qc_pass <- subset(hcc_11_qc, subset= hcc_11_qc$最终是否可用 ==1)
hcc_28_qc_pass <- subset(hcc_28_qc, subset= hcc_28_qc$最终是否可用 ==1)
hcc_29_qc_pass <- subset(hcc_29_qc, subset= hcc_29_qc$最终是否可用 ==1)


hcc3_base <- data.frame(var1="",var2="",var3="")[-1,]
hcc4_base <- data.frame(var1="",var2="",var3="")[-1,]
hcc7_base <- data.frame(var1="",var2="",var3="")[-1,]
hcc11_base <- data.frame(var1="",var2="",var3="")[-1,]
hcc28_base <- data.frame(var1="",var2="",var3="")[-1,]
hcc29_base <- data.frame(var1="",var2="",var3="")[-1,]




sample_list <- as.character(c("hcc3","hcc4","hcc7","hcc11","hcc28","hcc29"))
qc_list <- list(hcc_3_qc_pass,hcc_4_qc_pass,hcc_7_qc_pass,hcc_11_qc_pass,hcc_28_qc_pass,hcc_29_qc_pass)

for(i in 1:length(sample_list)){
  datapath <- paste("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/trio_seq/scMethy/methy_data/methy_file/",sample_list[i],"/CpG_100kb/",sep = "")
  setwd(datapath)
  getwd()
  fs <- list.files('./',pattern = "*.tsv")
  fs_pass <- as.character(1:nrow(qc_list[[i]]))
  j <- 1
  for (n in 1:length(fs)) {
    id <- str_split(paste(fs[n]),".CpG")[[1]][1]
    if(id %in% qc_list[[i]]$cell){
      fs_pass[j] <- fs[n]
      j <- j+1
    }
  }
  
  for(n in 1:length(fs_pass)){
    data <- fread(fs_pass[n])[,1:2]
    id <- str_split(paste(fs_pass[n]),".CpG")[[1]][1]
    colnames(data) <- c("pos",paste(id))
    if(n == 1){
      base <- data
      base <- as.data.frame(base)
    }
    else{
      base <- merge(base,data,by.x = 'pos',by.y = 'pos',all=TRUE, sort=TRUE)
    }
  }
  row.names(base) <- base[,1]
  base <- base[,-1]
  setwd("~/projects/hcc/analysis/trio_seq/scMethy")
  write.table(base,paste(sample_list[i],".txt",sep = ""))
}



##########hcc3#########
hcc3_base <- as.data.frame(fread("hcc3.txt"))
row.names(hcc3_base) <- hcc3_base$V1
hcc3_base <- hcc3_base[,-1]
hcc3_base_removena <- na.omit(hcc3_base)

hcc3_colanno <- as.data.frame(colnames(hcc3_base_removena))
hcc3_colanno[,2:3] <- str_split_fixed(hcc3_colanno[,1],"_",2)
colnames(hcc3_colanno) <- c("cellid","sample")
row.names(hcc3_colanno) <- hcc3_colanno[,1]
hcc3_colanno <- hcc3_colanno[,-3]
hcc3_colanno <- select(hcc3_colanno,sample)


hcc3_graph <- as.data.frame(str_split_fixed(row.names(hcc3_base_removena),"_",2))
colnames(hcc3_graph) <- c("chr","start")
hcc3_graph$pos <- row.names(hcc3_base_removena)
hcc3_graph$start <- as.numeric(hcc3_graph$start)
hcc3_graph$end <- hcc3_graph$start+1
hcc3_graph$start <- hcc3_graph$start*100000
hcc3_graph$end <- hcc3_graph$end*100000
hcc3_graph$end <- as.numeric(hcc3_graph$end)
hcc3_graph$start <- as.numeric(hcc3_graph$start)
hcc3_graph <- select(hcc3_graph,chr,start,end,pos)
write.table(hcc3_graph, "hcc3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

pwd_anno <- fread("big_pmd.bedGraph")
pwd_anno  <- select(pwd_anno ,V4,V8,V9)
pwd_anno <- aggregate(pwd_anno,by=list(pwd_anno$V4),max)
row.names_pwd_anno <- pwd_anno$Group.1 
pwd_anno <- select(pwd_anno,V4,V9)
pmd_anno <- as.data.frame(pwd_anno[,-1])
row.names(pmd_anno) <- row.names_pwd_anno
colnames(pmd_anno) <- "pmd_type"

hcc3_anno <- fread("hcc3_pmd.bedGraph")
hcc3_anno  <- select(hcc3_anno ,V4,V8,V9)
hcc3_anno <- aggregate(hcc3_anno,by=list(hcc3_anno$V4),max)
row.names_hcc3_anno <- hcc3_anno$Group.1 
hcc3_anno <- select(hcc3_anno,V4,V9)
hcc3_anno <- as.data.frame(hcc3_anno[,-1])
row.names(hcc3_anno) <- row.names_hcc3_anno
colnames(hcc3_anno) <- "meth_type"


hcc3_ann_colors=list(
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  sample=c('nt'='#DDE9C6','pt1'='#00ADC4','pt2'='#647BA8','pt3'='#4C1A72')
)

pheatmap(hcc3_base_removena,
         annotation_col = hcc3_colanno,annotation_row = hcc3_anno,
         annotation_colors = hcc3_ann_colors,
         show_rownames = F,show_colnames = F,         
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15,
         angle_col = 45)









##########hcc4##############
hcc4_base <- as.data.frame(fread("hcc4.txt"))
row.names(hcc4_base) <- hcc4_base$V1
hcc4_base <- hcc4_base[,-1]
hcc4_base_removena <- na.omit(hcc4_base)

hcc4_colanno <- as.data.frame(colnames(hcc4_base_removena))
hcc4_colanno[,2:3] <- str_split_fixed(hcc4_colanno[,1],"_",2)
colnames(hcc4_colanno) <- c("cellid","sample")
row.names(hcc4_colanno) <- hcc4_colanno[,1]
hcc4_colanno <- hcc4_colanno[,-3]
hcc4_colanno <- select(hcc4_colanno,sample)


hcc4_graph <- as.data.frame(str_split_fixed(row.names(hcc4_base_removena),"_",2))
colnames(hcc4_graph) <- c("chr","start")
hcc4_graph$pos <- row.names(hcc4_base_removena)
hcc4_graph$start <- as.numeric(hcc4_graph$start)
hcc4_graph$end <- hcc4_graph$start+1
hcc4_graph$start <- hcc4_graph$start*100000
hcc4_graph$end <- hcc4_graph$end*100000
hcc4_graph$end <- as.numeric(hcc4_graph$end)
hcc4_graph$start <- as.numeric(hcc4_graph$start)
hcc4_graph <- select(hcc4_graph,chr,start,end,pos)
write.table(hcc4_graph, "hcc4.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



hcc4_rowanno <- fread("hcc4_time.bedGraph")
hcc4_rowanno <- select(hcc4_rowanno,V4,V8,V9)
hcc4_rowanno <- aggregate(hcc4_rowanno,by=list(hcc4_rowanno$V8),max)
row.names(hcc4_rowanno) <- hcc4_rowanno$Group.1 
hcc4_rowanno <- select(hcc4_rowanno,V4)
colnames(hcc4_rowanno) <- "rep_time"

hcc4_anno <- fread("hcc4_pmd.bedGraph")
hcc4_anno  <- select(hcc4_anno ,V4,V8,V9)
hcc4_anno <- aggregate(hcc4_anno,by=list(hcc4_anno$V4),max)
row.names_hcc4_anno <- hcc4_anno$Group.1 
hcc4_anno <- select(hcc4_anno,V4,V9)
hcc4_anno <- as.data.frame(hcc4_anno[,-1])
row.names(hcc4_anno) <- row.names_hcc4_anno
colnames(hcc4_anno) <- "meth_type"

hcc4_ann_colors=list(
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  sample=c('nt'='#DDE9C6','pt1'='#00ADC4','pt2'='#647BA8','pt3'='#4C1A72')
)

pheatmap(hcc4_base_removena,
         annotation_col = hcc4_colanno,annotation_row = hcc4_anno,
         annotation_colors = hcc4_ann_colors,
         show_rownames = F,show_colnames = F,         
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15,
         angle_col = 45)









##########hcc7##############
hcc7_base <- as.data.frame(fread("hcc7.txt"))
row.names(hcc7_base) <- hcc7_base$V1
hcc7_base <- hcc7_base[,-1]
hcc7_base_removena <- na.omit(hcc7_base)

hcc7_colanno <- as.data.frame(colnames(hcc7_base_removena))
hcc7_colanno[,2:3] <- str_split_fixed(hcc7_colanno[,1],"_",2)
colnames(hcc7_colanno) <- c("cellid","sample")
row.names(hcc7_colanno) <- hcc7_colanno[,1]
hcc7_colanno <- hcc7_colanno[,-3]
hcc7_colanno <- select(hcc7_colanno,sample)


hcc7_graph <- as.data.frame(str_split_fixed(row.names(hcc7_base_removena),"_",2))
colnames(hcc7_graph) <- c("chr","start")
hcc7_graph$pos <- row.names(hcc7_base_removena)
hcc7_graph$start <- as.numeric(hcc7_graph$start)
hcc7_graph$end <- hcc7_graph$start+1
hcc7_graph$start <- hcc7_graph$start*100000
hcc7_graph$end <- hcc7_graph$end*100000
hcc7_graph$end <- as.numeric(hcc7_graph$end)
hcc7_graph$start <- as.numeric(hcc7_graph$start)
hcc7_graph <- select(hcc7_graph,chr,start,end,pos)
write.table(hcc7_graph, "hcc7.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



hcc7_rowanno <- fread("hcc7_time.bedGraph")
hcc7_rowanno <- select(hcc7_rowanno,V4,V8,V9)
hcc7_rowanno <- aggregate(hcc7_rowanno,by=list(hcc7_rowanno$V8),max)
row.names(hcc7_rowanno) <- hcc7_rowanno$Group.1 
hcc7_rowanno <- select(hcc7_rowanno,V4)
colnames(hcc7_rowanno) <- "rep_time"

hcc7_anno <- fread("hcc7_pmd.bedGraph")
hcc7_anno  <- select(hcc7_anno ,V4,V8,V9)
hcc7_anno <- aggregate(hcc7_anno,by=list(hcc7_anno$V4),max)
row.names_hcc7_anno <- hcc7_anno$Group.1 
hcc7_anno <- select(hcc7_anno,V4,V9)
hcc7_anno <- as.data.frame(hcc7_anno[,-1])
row.names(hcc7_anno) <- row.names_hcc7_anno
colnames(hcc7_anno) <- "meth_type"

hcc7_ann_colors=list(
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  sample=c('nt'='#DDE9C6','pt1'='#00ADC4','pt2'='#647BA8','pt4'='#4C1A72','pt5'='#A099C9')
)

pheatmap(hcc7_base_removena,
         annotation_col = hcc7_colanno,annotation_row = hcc7_anno,
         annotation_colors = hcc7_ann_colors,
         show_rownames = F,show_colnames = F,         
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15,
         angle_col = 45)




##########hcc11##############
hcc11_base <- as.data.frame(fread("hcc11.txt"))
row.names(hcc11_base) <- hcc11_base$V1
hcc11_base <- hcc11_base[,-1]
hcc11_base_removena <- na.omit(hcc11_base)

hcc11_colanno <- as.data.frame(colnames(hcc11_base_removena))
hcc11_colanno[,2:3] <- str_split_fixed(hcc11_colanno[,1],"_",2)
colnames(hcc11_colanno) <- c("cellid","sample")
row.names(hcc11_colanno) <- hcc11_colanno[,1]
hcc11_colanno <- hcc11_colanno[,-3]
hcc11_colanno <- select(hcc11_colanno,sample)


hcc11_graph <- as.data.frame(str_split_fixed(row.names(hcc11_base_removena),"_",2))
colnames(hcc11_graph) <- c("chr","start")
hcc11_graph$pos <- row.names(hcc11_base_removena)
hcc11_graph$start <- as.numeric(hcc11_graph$start)
hcc11_graph$end <- hcc11_graph$start+1
hcc11_graph$start <- hcc11_graph$start*100000
hcc11_graph$end <- hcc11_graph$end*100000
hcc11_graph$end <- as.numeric(hcc11_graph$end)
hcc11_graph$start <- as.numeric(hcc11_graph$start)
hcc11_graph <- select(hcc11_graph,chr,start,end,pos)
write.table(hcc11_graph, "hcc11.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



hcc11_rowanno <- fread("hcc11_time.bedGraph")
hcc11_rowanno <- select(hcc11_rowanno,V4,V8,V9)
hcc11_rowanno <- aggregate(hcc11_rowanno,by=list(hcc11_rowanno$V8),max)
row.names(hcc11_rowanno) <- hcc11_rowanno$Group.1 
hcc11_rowanno <- select(hcc11_rowanno,V4)
colnames(hcc11_rowanno) <- "rep_time"

hcc11_anno <- fread("hcc11_pmd.bedGraph")
hcc11_anno  <- select(hcc11_anno ,V4,V8,V9)
hcc11_anno <- aggregate(hcc11_anno,by=list(hcc11_anno$V4),max)
row.names_hcc11_anno <- hcc11_anno$Group.1 
hcc11_anno <- select(hcc11_anno,V4,V9)
hcc11_anno <- as.data.frame(hcc11_anno[,-1])
row.names(hcc11_anno) <- row.names_hcc11_anno
colnames(hcc11_anno) <- "meth_type"



hcc11_ann_colors=list(
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  sample=c('nt'='#DDE9C6','pt1'='#00ADC4','pt2'='#647BA8','pt3'='#A860A8','pt4'='#4C1A72','pt5'='#A099C9','pt6'='#AEE4F0')
)

pheatmap(hcc11_base_removena,
         annotation_col = hcc11_colanno,annotation_row = hcc11_anno,
         annotation_colors = hcc11_ann_colors,
         show_rownames = F,show_colnames = F,         
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15,
         angle_col = 45)







##########hcc28#####
hcc28_base <- as.data.frame(fread("hcc28.txt"))
row.names(hcc28_base) <- hcc28_base$V1
hcc28_base <- hcc28_base[,-1]
hcc28_base_removena <- na.omit(hcc28_base)

hcc28_colanno <- as.data.frame(colnames(hcc28_base_removena))
hcc28_colanno[,2:3] <- str_split_fixed(hcc28_colanno[,1],"_",2)
colnames(hcc28_colanno) <- c("cellid","sample")
row.names(hcc28_colanno) <- hcc28_colanno[,1]
hcc28_colanno <- hcc28_colanno[,-3]
hcc28_colanno <- select(hcc28_colanno,sample)


hcc28_graph <- as.data.frame(str_split_fixed(row.names(hcc28_base_removena),"_",2))
colnames(hcc28_graph) <- c("chr","start")
hcc28_graph$pos <- row.names(hcc28_base_removena)
hcc28_graph$start <- as.numeric(hcc28_graph$start)
hcc28_graph$end <- hcc28_graph$start+1
hcc28_graph$start <- hcc28_graph$start*100000
hcc28_graph$end <- hcc28_graph$end*100000
hcc28_graph$end <- as.numeric(hcc28_graph$end)
hcc28_graph$start <- as.numeric(hcc28_graph$start)
hcc28_graph <- select(hcc28_graph,chr,start,end,pos)
write.table(hcc28_graph, "hcc28.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



hcc28_rowanno <- fread("hcc28_time.bedGraph")
hcc28_rowanno <- select(hcc28_rowanno,V4,V8,V9)
hcc28_rowanno <- aggregate(hcc28_rowanno,by=list(hcc28_rowanno$V8),max)
row.names(hcc28_rowanno) <- hcc28_rowanno$Group.1 
hcc28_rowanno <- select(hcc28_rowanno,V4)
colnames(hcc28_rowanno) <- "rep_time"

hcc28_anno <- fread("hcc28_pmd.bedGraph")
hcc28_anno  <- select(hcc28_anno ,V4,V8,V9)
hcc28_anno <- aggregate(hcc28_anno,by=list(hcc28_anno$V4),max)
row.names_hcc28_anno <- hcc28_anno$Group.1 
hcc28_anno <- select(hcc28_anno,V4,V9)
hcc28_anno <- as.data.frame(hcc28_anno[,-1])
row.names(hcc28_anno) <- row.names_hcc28_anno
colnames(hcc28_anno) <- "meth_type"

hcc28_ann_colors=list(
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  sample=c('nt'='#DDE9C6','pt1'='#00ADC4','pt2'='#647BA8','pt4'='#4C1A72')
)

pheatmap(hcc28_base_removena,
         annotation_col = hcc28_colanno,annotation_row = hcc28_anno,
         annotation_colors = hcc28_ann_colors,
         show_rownames = F,show_colnames = F,         
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15,
         angle_col = 45)
















##########hcc29#####
hcc29_base <- as.data.frame(fread("hcc29.txt"))
row.names(hcc29_base) <- hcc29_base$V1
hcc29_base <- hcc29_base[,-1]
hcc29_base_removena <- na.omit(hcc29_base)



hcc29_colanno <- as.data.frame(colnames(hcc29_base_removena))
hcc29_colanno[,2:3] <- str_split_fixed(hcc29_colanno[,1],"_",2)
colnames(hcc29_colanno) <- c("cellid","sample")
row.names(hcc29_colanno) <- hcc29_colanno[,1]
hcc29_colanno <- hcc29_colanno[,-3]
hcc29_colanno <- select(hcc29_colanno,sample)



hcc29_graph <- as.data.frame(str_split_fixed(row.names(hcc29_base_removena),"_",2))
colnames(hcc29_graph) <- c("chr","start")
hcc29_graph$pos <- row.names(hcc29_base_removena)
hcc29_graph$start <- as.numeric(hcc29_graph$start)
hcc29_graph$end <- hcc29_graph$start+1
hcc29_graph$start <- hcc29_graph$start*100000
hcc29_graph$end <- hcc29_graph$end*100000
hcc29_graph$end <- as.numeric(hcc29_graph$end)
hcc29_graph$start <- as.numeric(hcc29_graph$start)
hcc29_graph <- select(hcc29_graph,chr,start,end,pos)
write.table(hcc29_graph, "hcc29.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


hcc29_rowanno <- fread("hcc29_time.bedGraph")
hcc29_rowanno <- select(hcc29_rowanno,V4,V8,V9)
hcc29_rowanno <- aggregate(hcc29_rowanno,by=list(hcc29_rowanno$V8),max)
row.names(hcc29_rowanno) <- hcc29_rowanno$Group.1 
hcc29_rowanno <- select(hcc29_rowanno,V4)
colnames(hcc29_rowanno) <- "rep_time"


hcc29_anno <- fread("hcc29_pmd.bedGraph")
hcc29_anno  <- select(hcc29_anno ,V4,V8,V9)
hcc29_anno <- aggregate(hcc29_anno,by=list(hcc29_anno$V4),max)
row.names_hcc29_anno <- hcc29_anno$Group.1 
hcc29_anno <- select(hcc29_anno,V4,V9)
hcc29_anno <- as.data.frame(hcc29_anno[,-1])
row.names(hcc29_anno) <- row.names_hcc29_anno
colnames(hcc29_anno) <- "meth_type"

hcc29_ann_colors=list(
  meth_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  sample=c('pt1'='#00ADC4','pt3'='#822994','pt4'='#4C1A72')
)

pheatmap(hcc29_base_removena,
         annotation_col = hcc29_colanno,annotation_row = hcc29_anno,
         annotation_colors = hcc29_ann_colors,
         show_rownames = F,show_colnames = F,
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15)






hcc3_for_merge <- hcc3_base_removena
colnames(hcc3_for_merge) <- paste("hcc3",colnames(hcc3_for_merge),sep = "_")

hcc4_for_merge <- hcc4_base_removena
colnames(hcc4_for_merge) <- paste("hcc4",colnames(hcc4_for_merge),sep = "_")

hcc7_for_merge <- hcc7_base_removena
colnames(hcc7_for_merge) <- paste("hcc7",colnames(hcc7_for_merge),sep = "_")

hcc11_for_merge <- hcc11_base_removena
colnames(hcc11_for_merge) <- paste("hcc11",colnames(hcc11_for_merge),sep = "_")

hcc28_for_merge <- hcc28_base_removena
colnames(hcc28_for_merge) <- paste("hcc28",colnames(hcc28_for_merge),sep = "_")

hcc29_for_merge <- hcc29_base_removena
colnames(hcc29_for_merge) <- paste("hcc29",colnames(hcc29_for_merge),sep = "_")



big_meth <- merge(hcc3_for_merge,hcc4_for_merge,by="row.names")
row.names(big_meth) <- big_meth$Row.names
big_meth <- big_meth[,-1]
big_meth <- merge(big_meth,hcc7_for_merge,by="row.names")
row.names(big_meth) <- big_meth$Row.names
big_meth <- big_meth[,-1]
big_meth <- merge(big_meth,hcc11_for_merge,by="row.names")
row.names(big_meth) <- big_meth$Row.names
big_meth <- big_meth[,-1]
big_meth <- merge(big_meth,hcc28_for_merge,by="row.names")
row.names(big_meth) <- big_meth$Row.names
big_meth <- big_meth[,-1]
big_meth <- merge(big_meth,hcc29_for_merge,by="row.names")
row.names(big_meth) <- big_meth$Row.names
big_meth <- big_meth[,-1]


hcc3_colanno_formerge <- hcc3_colanno
hcc3_colanno_formerge$sample <- paste("hcc3",hcc3_colanno_formerge$sample,sep = "_")
hcc3_colanno_formerge$origin <- "hcc3"
rownames(hcc3_colanno_formerge) <- paste("hcc3",rownames(hcc3_colanno_formerge),sep = "_")

hcc4_colanno_formerge <- hcc4_colanno
hcc4_colanno_formerge$sample <- paste("hcc4",hcc4_colanno_formerge$sample,sep = "_")
hcc4_colanno_formerge$origin <- "hcc4"
rownames(hcc4_colanno_formerge) <- paste("hcc4",rownames(hcc4_colanno_formerge),sep = "_")

hcc7_colanno_formerge <- hcc7_colanno
hcc7_colanno_formerge$sample <- paste("hcc7",hcc7_colanno_formerge$sample,sep = "_")
hcc7_colanno_formerge$origin <- "hcc7"
rownames(hcc7_colanno_formerge) <- paste("hcc7",rownames(hcc7_colanno_formerge),sep = "_")

hcc11_colanno_formerge <- hcc11_colanno
hcc11_colanno_formerge$sample <- paste("hcc11",hcc11_colanno_formerge$sample,sep = "_")
hcc11_colanno_formerge$origin <- "hcc11"
rownames(hcc11_colanno_formerge) <- paste("hcc11",rownames(hcc11_colanno_formerge),sep = "_")

hcc28_colanno_formerge <- hcc28_colanno
hcc28_colanno_formerge$sample <- paste("hcc28",hcc28_colanno_formerge$sample,sep = "_")
hcc28_colanno_formerge$origin <- "hcc28"
rownames(hcc28_colanno_formerge) <- paste("hcc28",rownames(hcc28_colanno_formerge),sep = "_")

hcc29_colanno_formerge <- hcc29_colanno
hcc29_colanno_formerge$sample <- paste("hcc29",hcc29_colanno_formerge$sample,sep = "_")
hcc29_colanno_formerge$origin <- "hcc29"
rownames(hcc29_colanno_formerge) <- paste("hcc29",rownames(hcc29_colanno_formerge),sep = "_")

big_colanno <- rbind(hcc3_colanno_formerge,hcc4_colanno_formerge,hcc7_colanno_formerge,hcc11_colanno_formerge,hcc28_colanno_formerge,hcc29_colanno_formerge)
big_colanno[,3:4] <- str_split_fixed(big_colanno[,1],"_",2)
for (i in 1:nrow(big_colanno)) {
  if(big_colanno[i,4]=='nt'){
    big_colanno[i,3] <- 'nt'
  }
   else{
     big_colanno[i,3] <- 'tumor'
   } 
}
big_colanno <- big_colanno[,-4]
colnames(big_colanno) <- c("sample","origin","tissue")


big_rowanno <- fread("hcc29_time.bedGraph")
big_rowanno <- select(big_rowanno,V4,V8,V9)
big_rowanno <- aggregate(big_rowanno,by=list(big_rowanno$V8),max)
row.names(big_rowanno) <- big_rowanno$Group.1 
big_rowanno <- select(big_rowanno,V4)
colnames(big_rowanno) <- "rep_time"



big_graph <- as.data.frame(str_split_fixed(row.names(big_meth),"_",2))
colnames(big_graph) <- c("chr","start")
big_graph$pos <- row.names(big_meth)
big_graph$start <- as.numeric(big_graph$start)
big_graph$end <- big_graph$start+1
big_graph$start <- big_graph$start*100000
big_graph$end <- big_graph$end*100000
big_graph$end <- as.numeric(big_graph$end)
big_graph$start <- as.numeric(big_graph$start)
big_graph <- select(big_graph,chr,start,end,pos)
write.table(big_graph, "big.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

big_rowanno <- fread("big_time.bedGraph")
big_rowanno <- select(big_rowanno,V4,V8,V9)
big_rowanno <- aggregate(big_rowanno,by=list(big_rowanno$V8),max)
row.names(big_rowanno) <- big_rowanno$Group.1 
big_rowanno <- select(big_rowanno,V4)
colnames(big_rowanno) <- "rep_time"



big_ann_colors=list(
  pmd_type=c('PMD'='#99CCFF','HMD'='#FF0000'),
  tissue=c('tumor'='#D42E00','nt'='#087FBF'),
  origin=c('hcc28'='#00ADC4','hcc29'='#4C1A72','hcc3'='#fb9795','hcc4'='#ffe9df','hcc7'='#7c5e8c','hcc11'='#b0e7ea')
)


pheatmap(big_meth,
         annotation_col = big_colanno,annotation_row = pmd_anno,
         annotation_colors = big_ann_colors,
         show_rownames = F,show_colnames = F,
         clustering_method = "average",
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15)




pheatmap(big_meth,
         annotation_col = big_colanno,annotation_row = big_rowanno,
         annotation_colors = big_ann_colors,
         show_rownames = F,show_colnames = F,
         clustering_method = "mcquitty",
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15)


bigmeth_pmd <- subset(pmd_anno, subset = pmd_anno$pmd_type == "PMD")

bigmeth_pmd_colmeans <- big_meth[row.names(bigmeth_pmd),]
bigmeth_pmd_colmeans <- as.data.frame(colMeans(bigmeth_pmd_colmeans))
bigmeth_pmd_colmeans[,"cell_id"] <- rownames(bigmeth_pmd_colmeans)
bigmeth_pmd_colmeans[,3:5] <- str_split_fixed(bigmeth_pmd_colmeans$cell_id,"_",3)
bigmeth_pmd_colmeans$V4 <- paste(bigmeth_pmd_colmeans$V3,bigmeth_pmd_colmeans$V4,sep = "_")
colnames(bigmeth_pmd_colmeans) <- c("mean_meth","cell_id","patient","sample")
bigmeth_pmd_colmeans <- bigmeth_pmd_colmeans[,-5]

bigmeth_pmd_colmeans_sample <- aggregate(bigmeth_pmd_colmeans$mean_meth,list(sample = bigmeth_pmd_colmeans$sample),mean)
bigmeth_pmd_colmeans_sample[,3:4] <- str_split_fixed(bigmeth_pmd_colmeans_sample$sample,"_",2)
bigmeth_pmd_colmeans_sample <- bigmeth_pmd_colmeans_sample[,-4]
colnames(bigmeth_pmd_colmeans_sample) <- c("sample","mean_meth","patient")


bigmeth_pmd_colmeans$sample = factor(bigmeth_pmd_colmeans$sample, levels=c('hcc3_nt','hcc3_pt1','hcc3_pt2','hcc3_pt3','hcc3_pt4',
                                                                           'hcc4_nt','hcc4_pt1','hcc4_pt2','hcc4_pt3',
                                                                           'hcc7_nt','hcc7_pt1','hcc7_pt2','hcc7_pt4','hcc7_pt5',
                                                                           'hcc11_nt','hcc11_pt1','hcc11_pt2','hcc11_pt3','hcc11_pt4',
                                                                           'hcc11_pt5','hcc11_pt6',
                                                                           'hcc28_nt','hcc28_pt1','hcc28_pt2','hcc28_pt4',
                                                                           'hcc29_nt','hcc29_pt1','hcc29_pt3','hcc29_pt4'))
bigmeth_pmd_colmeans$sample = factor(bigmeth_pmd_colmeans$sample, levels=sample_level)


bigmeth_pmd_colmeans[,5:6] <- str_split_fixed(bigmeth_pmd_colmeans$sample,"_",2)
for(i in 1:nrow(bigmeth_pmd_colmeans)){
  if(bigmeth_pmd_colmeans[i,6] == 'nt'){
    bigmeth_pmd_colmeans[i,6] <- 'normal'
  }
  else{
    bigmeth_pmd_colmeans[i,6] <- 'tumor'
  }
}

colnames(bigmeth_pmd_colmeans) <- c("mean_meth","cell_id","patient","sample","origin","label")


bigmeth_colmeans <- as.data.frame(colMeans(big_meth))
bigmeth_colmeans[,"cell_id"] <- rownames(bigmeth_colmeans)
bigmeth_colmeans[,3:5] <- str_split_fixed(bigmeth_colmeans$cell_id,"_",3)
bigmeth_colmeans$V4 <- paste(bigmeth_colmeans$V3,bigmeth_colmeans$V4,sep = "_")
colnames(bigmeth_colmeans) <- c("mean_meth","cell_id","patient","sample")
bigmeth_colmeans <- bigmeth_colmeans[,-5]

bigmeth_colmeans_sample <- aggregate(bigmeth_colmeans$mean_meth,list(sample = bigmeth_colmeans$sample),mean)
bigmeth_colmeans_sample[,3:4] <- str_split_fixed(bigmeth_colmeans_sample$sample,"_",2)
bigmeth_colmeans_sample <- bigmeth_colmeans_sample[,-4]
colnames(bigmeth_colmeans_sample) <- c("sample","mean_meth","patient")


bigmeth_colmeans$sample = factor(bigmeth_colmeans$sample, levels=c('hcc3_nt','hcc3_pt1','hcc3_pt2','hcc3_pt3','hcc3_pt4',
                                                                   'hcc4_nt','hcc4_pt1','hcc4_pt2','hcc4_pt3',
                                                                   'hcc7_nt','hcc7_pt1','hcc7_pt2','hcc7_pt4','hcc7_pt5',
                                                                   'hcc11_nt','hcc11_pt1','hcc11_pt2','hcc11_pt3','hcc11_pt4',
                                                                   'hcc11_pt5','hcc11_pt6',
                                                                   'hcc28_nt','hcc28_pt1','hcc28_pt2','hcc28_pt4',
                                                                   'hcc29_nt','hcc29_pt1','hcc29_pt3','hcc29_pt4'))
sample_level <- rev(c('hcc3_nt','hcc3_pt1','hcc3_pt2','hcc3_pt3','hcc3_pt4',
                  'hcc4_nt','hcc4_pt1','hcc4_pt2','hcc4_pt3',
                  'hcc7_nt','hcc7_pt1','hcc7_pt2','hcc7_pt4','hcc7_pt5',
                  'hcc11_nt','hcc11_pt1','hcc11_pt2','hcc11_pt3','hcc11_pt4',
                  'hcc11_pt5','hcc11_pt6',
                  'hcc28_nt','hcc28_pt1','hcc28_pt2','hcc28_pt4',
                  'hcc29_nt','hcc29_pt1','hcc29_pt3','hcc29_pt4'))
bigmeth_colmeans$sample = factor(bigmeth_colmeans$sample, levels=sample_level)


bigmeth_colmeans[,5:6] <- str_split_fixed(bigmeth_colmeans$sample,"_",2)
for(i in 1:nrow(bigmeth_colmeans)){
  if(bigmeth_colmeans[i,6] == 'nt'){
    bigmeth_colmeans[i,6] <- 'normal'
  }
  else{
    bigmeth_colmeans[i,6] <- 'tumor'
  }
}

colnames(bigmeth_colmeans) <- c("mean_meth","cell_id","patient","sample","origin","label")

ggplot(bigmeth_colmeans,aes(x=patient, y=mean_meth,fill=origin))+
  geom_flat_violin(scale = "width",trim = F)+coord_flip()+geom_jitter(width = 0.1,size=0.5)
ggplot(bigmeth_pmd_colmeans,aes(x=patient, y=mean_meth,fill=label))+
  geom_flat_violin(scale = "width",trim = F)+coord_flip()+geom_jitter(width = 0.1,size=0.5)

ggplot(bigmeth_colmeans,aes(x=sample, y=mean_meth,fill=origin))+
  geom_flat_violin(scale = "width",trim = F)+geom_jitter(width = 0.1,size=0.5)+coord_flip()
ggplot(bigmeth_pmd_colmeans,aes(x=sample, y=mean_meth,fill=label))+
  geom_flat_violin(scale = "width",trim = F)+geom_jitter(width = 0.1,size=0.5)+coord_flip()

ggplot(bigmeth_colmeans_sample,aes(x=patient, y=mean_meth,fill=patient))+
  geom_flat_violin(trim = F)+NoLegend()+coord_flip()

ggplot(bigmeth_colmeans,aes(x=patient, y=mean_meth,fill=sample))+
  geom_flat_violin(trim = F)+NoLegend()+coord_flip()

ggplot(bigmeth_colmeans,aes(x=patient, y=mean_meth,fill=patient))+
  geom_flat_violin(trim = F)+NoLegend()+coord_flip()+

ggplot(bigmeth_colmeans,aes(x=patient, y=mean_meth,fill=sample))+
  geom_boxplot()+NoLegend()+coord_flip()

bigmeth_colmeans_hcc11 <- subset(bigmeth_colmeans,subset = patient=="hcc11")
ggplot(bigmeth_colmeans_hcc11,aes(x=sample, y=mean_meth,fill=sample))+
  geom_flat_violin(scale = "width")+coord_flip()+geom_jitter(width = 0.1,size=0.5)


bigmeth_colmeans_hcc3 <- subset(bigmeth_colmeans,subset = patient=="hcc3")
ggplot(bigmeth_colmeans_hcc3,aes(x=sample, y=mean_meth,fill=sample))+
  geom_flat_violin(scale = "width")+coord_flip()+geom_jitter(width = 0.1,size=0.5)



bigmeth_colmeans_hcc4 <- subset(bigmeth_colmeans,subset = patient=="hcc4")
ggplot(bigmeth_colmeans_hcc4,aes(x=sample, y=mean_meth,fill=sample))+
  geom_flat_violin(scale = "width")+coord_flip()+geom_jitter(width = 0.1,size=0.5)


bigmeth_colmeans_hcc7 <- subset(bigmeth_colmeans,subset = patient=="hcc7")
ggplot(bigmeth_colmeans_hcc7,aes(x=sample, y=mean_meth,fill=sample))+
  geom_flat_violin(scale = "width")+coord_flip()+geom_jitter(width = 0.1,size=0.5)


bigmeth_colmeans_hcc28 <- subset(bigmeth_colmeans,subset = patient=="hcc28")
ggplot(bigmeth_colmeans_hcc28,aes(x=sample, y=mean_meth,fill=sample))+
  geom_flat_violin(scale = "width")+coord_flip()+geom_jitter(width = 0.1,size=0.5)


bigmeth_colmeans_hcc29 <- subset(bigmeth_colmeans,subset = patient=="hcc29")
ggplot(bigmeth_colmeans_hcc29,aes(x=sample, y=mean_meth,fill=sample))+
  geom_flat_violin(scale = "width")+coord_flip()+geom_jitter(width = 0.1,size=0.5)




