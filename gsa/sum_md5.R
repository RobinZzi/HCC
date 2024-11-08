 
library(tools)

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
  
  cbind(general_df,r2_df)
}  


result <- listFilesInDirectory()  

