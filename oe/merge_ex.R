


setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge")


setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex1/meth_ex/methy_report")
EX1_G6_OE1 <- fread("HepG2-SNHG6-OE1_remove_dup.CpG_report.txt")
EX1_G6_OE2 <- fread("HepG2-SNHG6-OE2_remove_dup.CpG_report.txt")
EX1_G6_OE3 <- fread("HepG2-SNHG6-OE3_remove_dup.CpG_report.txt")
EX1_G6_C1 <- fread("HepG2-SNHG6-C1_remove_dup.CpG_report.txt")
EX1_G6_C2 <- fread("HepG2-SNHG6-C2_remove_dup.CpG_report.txt")
EX1_G6_C3 <- fread("HepG2-SNHG6-C3_remove_dup.CpG_report.txt")


EX1_GA_OE1 <- fread("HepG2-GADD45A-OE1_remove_dup.CpG_report.txt")
EX1_GA_OE2 <- fread("HepG2-GADD45A-OE2_remove_dup.CpG_report.txt")
EX1_GA_OE3 <- fread("HepG2-GADD45A-OE3_remove_dup.CpG_report.txt")
EX1_GA_C1 <- fread("HepG2-GADD45A-C1_remove_dup.CpG_report.txt")
EX1_GA_C2 <- fread("HepG2-GADD45A-C2_remove_dup.CpG_report.txt")
EX1_GA_C3 <- fread("HepG2-GADD45A-C3_remove_dup.CpG_report.txt")


setwd("~/projects/hcc/analysis/oe/ex2/methy_merge/report")
EX2_G6_OE1 <- fread("SNHG6-OE1.deduplicated.CpG_report.txt")
EX2_G6_OE2 <- fread("SNHG6-OE2.deduplicated.CpG_report.txt")
EX2_G6_OE3 <- fread("SNHG6-OE3.deduplicated.CpG_report.txt")
EX2_G6_C1 <- fread("SNHG6-C1.deduplicated.CpG_report.txt")
EX2_G6_C2 <- fread("SNHG6-C2.deduplicated.CpG_report.txt")
EX2_G6_C3 <- fread("SNHG6-C3.deduplicated.CpG_report.txt")


EX2_GA_OE1 <- fread("GADD45A-OE1.deduplicated.CpG_report.txt")
EX2_GA_OE2 <- fread("GADD45A-OE2.deduplicated.CpG_report.txt")
EX2_GA_OE3 <- fread("GADD45A-OE3.deduplicated.CpG_report.txt")
EX2_GA_C1 <- fread("GADD45A-C1.deduplicated.CpG_report.txt")
EX2_GA_C2 <- fread("GADD45A-C2.deduplicated.CpG_report.txt")
EX2_GA_C3 <- fread("GADD45A-C3.deduplicated.CpG_report.txt")


setwd("~/projects/hcc/analysis/oe/merge")
G6_OE1_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(G6_OE1_SUM) <- c("chr","cite","meth","demeth")
G6_OE1_SUM$meth <- EX1_G6_OE1$V4 + EX2_G6_OE1$V4
G6_OE1_SUM$demeth <- EX1_G6_OE1$V5 + EX2_G6_OE1$V5
G6_OE1_SUM[,1:2] <- EX1_G6_OE1[,1:2]
write.table(G6_OE1_SUM,"G6_OE1_SUM.txt")



G6_OE2_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(G6_OE2_SUM) <- c("chr","cite","meth","demeth")
G6_OE2_SUM$meth <- EX1_G6_OE2$V4 + EX2_G6_OE2$V4
G6_OE2_SUM$demeth <- EX1_G6_OE2$V5 + EX2_G6_OE2$V5
G6_OE2_SUM[,1:2] <- EX1_G6_OE2[,1:2]
write.table(G6_OE2_SUM,"G6_OE2_SUM.txt")



G6_OE3_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(G6_OE3_SUM) <- c("chr","cite","meth","demeth")
G6_OE3_SUM$meth <- EX1_G6_OE3$V4 + EX2_G6_OE3$V4
G6_OE3_SUM$demeth <- EX1_G6_OE3$V5 + EX2_G6_OE3$V5
G6_OE3_SUM[,1:2] <- EX1_G6_OE3[,1:2]
write.table(G6_OE3_SUM,"G6_OE3_SUM.txt")



G6_C1_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(G6_C1_SUM) <- c("chr","cite","meth","demeth")
G6_C1_SUM$meth <- EX1_G6_C1$V4 + EX2_G6_C1$V4
G6_C1_SUM$demeth <- EX1_G6_C1$V5 + EX2_G6_C1$V5
G6_C1_SUM[,1:2] <- EX1_G6_C1[,1:2]
write.table(G6_C1_SUM,"G6_C1_SUM.txt")



G6_C2_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(G6_C2_SUM) <- c("chr","cite","meth","demeth")
G6_C2_SUM$meth <- EX1_G6_C2$V4 + EX2_G6_C2$V4
G6_C2_SUM$demeth <- EX1_G6_C2$V5 + EX2_G6_C2$V5
G6_C2_SUM[,1:2] <- EX1_G6_C2[,1:2]
write.table(G6_C2_SUM,"G6_C2_SUM.txt")



G6_C3_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(G6_C3_SUM) <- c("chr","cite","meth","demeth")
G6_C3_SUM$meth <- EX1_G6_C3$V4 + EX2_G6_C3$V4
G6_C3_SUM$demeth <- EX1_G6_C3$V5 + EX2_G6_C3$V5
G6_C3_SUM[,1:2] <- EX1_G6_C3[,1:2]
write.table(G6_C3_SUM,"G6_C3_SUM.txt")



GA_OE1_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(GA_OE1_SUM) <- c("chr","cite","meth","demeth")
GA_OE1_SUM$meth <- EX1_GA_OE1$V4 + EX2_GA_OE1$V4
GA_OE1_SUM$demeth <- EX1_GA_OE1$V5 + EX2_GA_OE1$V5
GA_OE1_SUM[,1:2] <- EX1_GA_OE1[,1:2]
write.table(GA_OE1_SUM,"GA_OE1_SUM.txt")



GA_OE2_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(GA_OE2_SUM) <- c("chr","cite","meth","demeth")
GA_OE2_SUM$meth <- EX1_GA_OE2$V4 + EX2_GA_OE2$V4
GA_OE2_SUM$demeth <- EX1_GA_OE2$V5 + EX2_GA_OE2$V5
GA_OE2_SUM[,1:2] <- EX1_GA_OE2[,1:2]
write.table(GA_OE2_SUM,"GA_OE2_SUM.txt")



GA_OE3_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(GA_OE3_SUM) <- c("chr","cite","meth","demeth")
GA_OE3_SUM$meth <- EX1_GA_OE3$V4 + EX2_GA_OE3$V4
GA_OE3_SUM$demeth <- EX1_GA_OE3$V5 + EX2_GA_OE3$V5
GA_OE3_SUM[,1:2] <- EX1_GA_OE3[,1:2]
write.table(GA_OE3_SUM,"GA_OE3_SUM.txt")



GA_C1_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(GA_C1_SUM) <- c("chr","cite","meth","demeth")
GA_C1_SUM$meth <- EX1_GA_C1$V4 + EX2_GA_C1$V4
GA_C1_SUM$demeth <- EX1_GA_C1$V5 + EX2_GA_C1$V5
GA_C1_SUM[,1:2] <- EX1_GA_C1[,1:2]
write.table(GA_C1_SUM,"GA_C1_SUM.txt")



GA_C2_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(GA_C2_SUM) <- c("chr","cite","meth","demeth")
GA_C2_SUM$meth <- EX1_GA_C2$V4 + EX2_GA_C2$V4
GA_C2_SUM$demeth <- EX1_GA_C2$V5 + EX2_GA_C2$V5
GA_C2_SUM[,1:2] <- EX1_GA_C2[,1:2]
write.table(GA_C2_SUM,"GA_C2_SUM.txt")



GA_C3_SUM <- as.data.frame(matrix(0,58607920,4))
colnames(GA_C3_SUM) <- c("chr","cite","meth","demeth")
GA_C3_SUM$meth <- EX1_GA_C3$V4 + EX2_GA_C3$V4
GA_C3_SUM$demeth <- EX1_GA_C3$V5 + EX2_GA_C3$V5
GA_C3_SUM[,1:2] <- EX1_GA_C3[,1:2]
write.table(GA_C3_SUM,"GA_C3_SUM.txt")

chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
          "chr16","chr17","chr18","chr19","chr20","chr21","chr22")

setwd("~/projects/hcc/analysis/oe/merge/merged_report")
fs <- list.files()
for (i in 1:length(fs)) {
  setwd("~/projects/hcc/analysis/oe/merge/merged_report")
  data <- fread(paste(fs[i]))
  data <- subset(data,data$chr %in% chrs)
  data <- data[,-1]
  data$sum <- round(data[,2]/100000)
  data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
  data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
  sample_id <-str_split_fixed(fs[i],".txt",2)[1]
  data[,paste(sample_id)] <- data[,2]/(data[,2]+data[,3])
  setwd("~/projects/hcc/analysis/oe/merge/merged_level")
  write.table(data,paste(sample_id,".txt",sep = ""))
}
 
save.image("merge_ex.Rdata")


