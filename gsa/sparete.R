bigseu <- readRDS("~/projects/hcc/analysis/scRNA_backup/bigseu.RDS")
bigseu <- readRDS("~/projects/hcc/analysis/merged_scrna/strt_merge.RDS")
table(bigseu$patient)
# HCC1  HCC2  HCC3  HCC4  HCC5  HCC6  HCC7  HCC8  HCC9 
#13761 21715  1778   948  1231  2081  3686 27767 25252

table(bigseu$lib.method)
#10x dropseq 
#88495    9724 
table(bigseu$sample_pt)

samplelist <- as.data.frame(table(bigseu$patient_pt))
samplelist <- as.data.frame(table(bigseu$sample_pt))
pt_list <- unlist(unlist(samplelist$Var1))

for (i in 1:length(pt_list)){
 setwd("~/projects/hcc/analysis/strt_back_to_10x")
 name <- pt_list[i]
 dir.create(paste(toupper(name)))
 setwd(paste("~/projects/hcc/analysis/strt_back_to_10x/",toupper(name),sep = ""))
 seu <-  subset(bigseu,subset = patient_pt == paste(name))
 ct=GetAssayData(object = seu, assay = "RNA", slot = "counts") 
 ct=as.data.frame(ct)
 write.table(data.frame(rownames(ct),rownames(ct)),file = 'genes.tsv',
            quote = F,sep = '\t',
            col.names = F,row.names = F)
 write.table(colnames(ct),file = 'barcodes.tsv',quote = F,
            col.names = F,row.names = F)
 file="matrix.mtx"
 sink(file)
 cat("%%MatrixMarket matrix coordinate integer general\n")
 cat("%\n")
 cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
 sink()
 tmp=ct[1:5,1:4]
 tmp
 tmp=do.call(rbind,lapply(1:ncol(ct),function(i){
  return(data.frame(row=1:nrow(ct),
                    col=i,
                    exp=ct[,i]))
  }) )
 tmp=tmp[tmp$exp>0,]
 head(tmp)
 write.table(tmp,file = 'matrix.mtx',quote = F,
            col.names = F,row.names = F,append = T )
 paste(name,"is done",sep = " ")
}
