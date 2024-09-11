#10-Mb genomic windows


setwd("~/projects/hcc/analysis/cell_re_rrbs/samtools_depth_result")




rm(list=ls())
folders <- list.files()
data_path <- "~/projects/hcc/analysis/cell_re_rrbs/samtools_depth_result"
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
          "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
for (i in 1:length(folders)) {
  setwd(paste(data_path,folders[i],sep = "/"))
  mtxs <- list.files()
  for(j in 1:length(mtxs)){
    sample_name <- str_split(mtxs[j],"_sorted.depth")[[1]][1]
    print(paste(sample_name,"start",sep = ":"))
    data <- fread(mtxs[j])
    data <- subset(data,subset=V1 %in% chrs)
    data <- mutate(data,bin_id=floor(V2/10000000+1))
    data <- mutate(data,bin_id=paste(V1,bin_id,sep = "_"))
    data <- aggregate(data$V3,by=list(data$bin_id),FUN=sum)
    print(paste(sample_name,"end",sep = ":"))
    write.table(data,paste(sample_name,"_10M",".txt",sep = ""))
  }
  
}




