rm(list=ls())
library(xlsx)
library(stringr)
library(pheatmap)
save.image("methy_merge.Rdata")
setwd("~/projects/hcc/data/trio_seq/scMethy")
hcc_28_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc28") 
hcc_29_qc <- read.xlsx2("methy.qc.xlsx",sheetName = "hcc29") 


hcc_28_qc_pass <- subset(hcc_28_qc, subset= hcc_28_qc$最终是否可用 ==1)
hcc_29_qc_pass <- subset(hcc_29_qc, subset= hcc_29_qc$最终是否可用 ==1)


setwd("~/projects/hcc/data/trio_seq/scMethy/methy_data/methy_file/hcc28/CpG_100kb")
fs_hcc28 <- list.files('./',pattern = "*.tsv")
fs_hcc28_pass <- as.character(1:124)
j <- 1
for (i in 1:length(fs_hcc28)) {
  id <- str_split(paste(fs_hcc28[i]),".CpG")[[1]][1]
  if(id %in% hcc_28_qc_pass$cell){
    fs_hcc28_pass[j] <- fs_hcc28[i]
    j <- j+1
  }
}

for(i in 1:length(fs_hcc28_pass)){
  data <- fread(fs_hcc28_pass[i])[,1:2]
  id <- str_split(paste(fs_hcc28_pass[i]),".CpG")[[1]][1]
  colnames(data) <- c("pos",paste(id))
  if(i == 1){
    hcc28_base <- data
    hcc28_base <- as.data.frame(hcc28_base)
  }
  else{
    hcc28_base <- merge(hcc28_base,data,by.x = 'pos',by.y = 'pos',all=TRUE, sort=TRUE)
  }

}
row.names(hcc28_base) <- hcc28_base[,1]
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



pheatmap(hcc28_base_removena,
         annotation_col = hcc28_colanno,annotation_row = hcc28_rowanno,
         annotation_colors = hcc28_ann_colors,
         show_rownames = F,show_colnames = F,         
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15,
         angle_col = 45)




hcc28_ann_colors=list(
  rep_time=c('LRD'='#0039B3','ERD'='#FF1962','DTZ'='#FFE0EC','UTZ'='#FFF599'),
  sample=c('nt'='#DDE9C6','pt1'='#00ADC4','pt2'='#647BA8','pt4'='#4C1A72')
  )

hcc29_ann_colors=list(
  rep_time=c('LRD'='#0039B3','ERD'='#FF1962','DTZ'='#FFE0EC','UTZ'='#FFF599'),
  sample=c('pt1'='#00ADC4','pt3'='#822994','pt4'='#4C1A72')
)










  setwd("~/projects/hcc/data/trio_seq/scMethy/methy_data/methy_file/hcc29/CpG_100kb")
  fs_hcc29 <- list.files('./',pattern = "*.tsv")
  fs_hcc29_pass <- as.character(1:nrow(hcc_29_qc_pass))
  j <- 1
  for (i in 1:length(fs_hcc29)) {
    id <- str_split(paste(fs_hcc29[i]),".CpG")[[1]][1]
    if(id %in% hcc_29_qc_pass$cell){
      fs_hcc29_pass[j] <- fs_hcc29[i]
      j <- j+1
    }
  }
  
  for(i in 1:length(fs_hcc29_pass)){
    data <- fread(fs_hcc29_pass[i])[,1:2]
    id <- str_split(paste(fs_hcc29_pass[i]),".CpG")[[1]][1]
    colnames(data) <- c("pos",paste(id))
    if(i == 1){
      hcc29_base <- data
      hcc29_base <- as.data.frame(hcc29_base)
    }
    else{
      hcc29_base <- merge(hcc29_base,data,by.x = 'pos',by.y = 'pos',all=TRUE, sort=TRUE)
    }
    
  }
  row.names(hcc29_base) <- hcc29_base[,1]
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
 
pheatmap(hcc29_base_removena,
         annotation_col = hcc29_colanno,annotation_row = hcc29_rowanno,
         annotation_colors = hcc29_ann_colors,
         show_rownames = F,show_colnames = F,
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15,
         angle_col = 180)


hcc29_for_merge <- hcc29_base_removena
colnames(hcc29_for_merge) <- paste("hcc29",colnames(hcc29_for_merge),sep = "_")

hcc28_for_merge <- hcc28_base_removena
colnames(hcc28_for_merge) <- paste("hcc28",colnames(hcc28_for_merge),sep = "_")

big_meth <- merge(hcc28_for_merge,hcc29_for_merge,by="row.names")
row.names(big_meth) <- big_meth$Row.names
big_meth <- big_meth[,-1]

hcc28_colanno_formerge <- hcc28_colanno
hcc28_colanno_formerge$sample <- paste("hcc28",hcc28_colanno_formerge$sample,sep = "_")
hcc28_colanno_formerge$origin <- "hcc28"
rownames(hcc28_colanno_formerge) <- paste("hcc28",rownames(hcc28_colanno_formerge),sep = "_")

hcc29_colanno_formerge <- hcc29_colanno
hcc29_colanno_formerge$sample <- paste("hcc29",hcc29_colanno_formerge$sample,sep = "_")
hcc29_colanno_formerge$origin <- "hcc29"
rownames(hcc29_colanno_formerge) <- paste("hcc29",rownames(hcc29_colanno_formerge),sep = "_")

big_colanno <- rbind(hcc29_colanno_formerge,hcc28_colanno_formerge)
big_rowanno <- subset(hcc28_rowanno,subset=row.names(hcc28_rowanno) %in% row.names(big_meth))



pheatmap(big_meth,
         annotation_col = big_colanno,annotation_row = big_rowanno,
         annotation_colors = big_ann_colors,
         show_rownames = F,show_colnames = F,
         clustering_method = "complete",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15)


big_ann_colors=list(
  rep_time=c('LRD'='#0039B3','ERD'='#FF1962','DTZ'='#FFE0EC','UTZ'='#FFF599'),
  origin=c('hcc28'='#00ADC4','hcc29'='#4C1A72')
)




