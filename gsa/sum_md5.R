 
library(tools)
library(stringr)
rm(list = ls())
listFilesInDirectory <- function() {  
  all_files <- list.files(pattern = "*.gz$", full.names = TRUE)  
  all_files <- sapply(all_files, function(x) sub("^\\./", "", x))

  general_list <- list()  
  r2_list <- list()  

  for (gz_file in all_files) {  
    md5_file <- paste0(gz_file,".md5")
    
    if (file.exists(md5_file)) {  
      md5_content <- readLines(md5_file, warn = FALSE)  
    } else {  
      md5_content <- NA  
    }  
    md5_content <- substr(md5_content, 1, 32)

    if (grepl("R2", gz_file)) {  
      r2_list[[gz_file]] <- md5_content  
    } else {  
      general_list[[gz_file]] <- md5_content  
    }  
  }  
  
 
  general_df <- data.frame(  
    gz_file = names(general_list),  
    md5_content = I(general_list),  
    stringsAsFactors = FALSE 
  )  
  
  r2_df <- data.frame(  
    r2_gz_file = names(r2_list),  
    r2_md5_content = I(r2_list),  
    stringsAsFactors = FALSE  
  )  
  
  list(general_df,r2_df)
}  
origin_dir <- "/mnt/transposon1/zhangyanxiaoLab/zhangliwen/hcc_multiomics_gsa_upload/"
origin_dir <- "/mnt/usb/"

data_dir <- paste(origin_dir,'data',sep='/')
md5_dir <- paste(origin_dir,'md5_table',sep='/')
setwd(data_dir)
folder_list <- list.files()
for (i in 1:length(folder_list)) {
  first_dir <- paste(data_dir,folder_list[i],sep = "/")
  setwd(first_dir)
  second_list <-list.files()
  for (j in 1:length(second_list)) {
    second_dir <- paste(first_dir,second_list[j],sep='/')
    setwd(second_dir)
  result <- listFilesInDirectory()  
    general_df <- result[[1]]  
    r2_df <- result[[2]] 
    if(nrow(r2_df) != 0){
      general_df <- cbind(general_df,r2_df)
    }
   if(folder_list[i] == 'scPBAT'){
    filenames <- rownames(general_df)
    t_numbers <- as.numeric(str_extract(filenames, "(?<=T)\\d+"))+1 
    t_numbers[is.na(t_numbers)] <- 1
    d_numbers <- as.numeric(str_extract(filenames, "(?<=D)\\d+"))
    sample_order <- t_numbers*1000+d_numbers
    general_df <- general_df[order(sample_order), ]  
    }
    setwd(md5_dir)
    write.csv(general_df,paste0(second_list[j],'_general.csv'))
  }
}
 

