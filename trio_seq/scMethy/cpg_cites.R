


library(data.table)
library(stringr)
rm(list=ls())
sample_list <- as.character(c("hcc3","hcc4","hcc7","hcc11","hcc28","hcc29"))

for(i in 1:length(sample_list)){
  datapath <- paste("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/trio_seq/scMethy/methy_data/methy_file/",sample_list[i],"/CpG_profile/",sep = "")
  setwd(datapath)
  fs_sample <- list.files('./')
  getwd()
  for (j in 1:length(fs_sample)) {
    datapath <- paste("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/trio_seq/scMethy/methy_data/methy_file/",sample_list[i],"/CpG_profile/",fs_sample[j],sep = "")
    setwd(datapath)
    fs <- list.files('./',pattern = "*.gz")
    for(n in 1:length(fs)){
      data <- fread(fs[n])
      id <- str_split(paste(fs[n]),".CpG")[[1]][1]
      colnames(data) <- c("chr","start","end","level","meth","demeth")
      data$pos <- paste(data$chr,data$start,sep = "_")
      data <- select(data,pos,meth,demeth)
      if(n == 1){
        base <- data
        base <- as.data.frame(base)
      }else{
        base <- merge(base,data,by.x = 'pos',by.y = 'pos',all=TRUE, sort=TRUE)
        base[is.na(base)] <- 0
        base$meth <- base$meth.x+base$meth.y
        base$demeth <- base$demeth.x+base$demeth.y
        base <- select(base,pos,meth,demeth)
      }
      print(id)
    }
    print(paste(sample_list[i],"_",fs_sample[j],".txt",sep = ""))
    setwd("~/projects/hcc/analysis/trio_seq/scMethy")
    write.table(base,paste(sample_list[i],"_",fs_sample[j],".txt",sep = ""))
  }
}
