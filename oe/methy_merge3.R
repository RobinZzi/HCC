setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report")

rm(list=ls())
library(data.table)
library(R.utils)
library(tidyr)
library(pheatmap)
library(edgeR)
library(stringr)
library(FactoMineR)
library(dplyr)


load("methy_merge.Rdata")
save.image("methy_merge.Rdata")

chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
          "chr16","chr17","chr18","chr19","chr20","chr21","chr22")
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report/report")
fs <- list.files()
for (i in 1:length(fs)) {
  setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report/report")
  data <- fread(paste(fs[i]))
  data <- subset(data,data$V1 %in% chrs)
  colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
  data[,3] <- round(data[,2]/100000)
  data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
  data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
  sample_id <-str_split_fixed(fs[i],"_report.txt",2)[1]
  data[,"count"] <- data[,2]+data[,3]
  data[,paste(sample_id)] <- data[,2]/(data[,2]+data[,3])
  setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report")
  write.table(data,paste(sample_id,".txt",sep = ""))
}


GA_C1 <- fread("GADD45A-C1.txt")
GA_C1[,"count"] <-GA_C1$meth +GA_C1$dmeth
GA_C2 <- fread("GADD45A-C2.txt")
GA_C2[,"count"] <-GA_C2$meth +GA_C2$dmeth
GA_C3 <- fread("GADD45A-C3.txt")
GA_C3[,"count"] <-GA_C3$meth +GA_C3$dmeth
GA_E1 <- fread("GADD45A-OE1.txt")
GA_E1[,"count"] <-GA_E1$meth +GA_E1$dmeth
GA_E2 <- fread("GADD45A-OE2.txt")
GA_E2[,"count"] <-GA_E2$meth +GA_E2$dmeth
GA_E3 <- fread("GADD45A-OE3.txt")
GA_E3[,"count"] <-GA_E3$meth +GA_E3$dmeth




GA_C1_fragment <- GA_C1$count[GA_C1$count <= 5000]
head(GA_C1_fragment)
GA_C1_fragment <- data.frame(GA_C1_fragment)
GA_C1_res <- hist(GA_C1_fragment$GA_C1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_C1_res$breaks),
     y = c(0, 0, GA_C1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_C1 Sample reads stat")

GA_C2_fragment <- GA_C2$count[GA_C2$count <= 5000]
head(GA_C2_fragment)
GA_C2_fragment <- data.frame(GA_C2_fragment)
GA_C2_res <- hist(GA_C2_fragment$GA_C2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_C2_res$breaks),
     y = c(0, 0, GA_C2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_C2 Sample reads stat")

GA_C3_fragment <- GA_C3$count[GA_C3$count <= 5000]
head(GA_C2_fragment)
GA_C3_fragment <- data.frame(GA_C3_fragment)
GA_C3_res <- hist(GA_C3_fragment$GA_C3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_C3_res$breaks),
     y = c(0, 0, GA_C3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_C3 Sample reads stat")

GA_E1_fragment <- GA_E1$count[GA_E1$count <= 5000]
head(GA_E1_fragment)
GA_E1_fragment <- data.frame(GA_E1_fragment)
GA_E1_res <- hist(GA_E1_fragment$GA_E1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_E1_res$breaks),
     y = c(0, 0, GA_E1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_E1 Sample reads stat")

GA_E2_fragment <- GA_E2$count[GA_E2$count <= 5000]
head(GA_E2_fragment)
GA_E2_fragment <- data.frame(GA_E2_fragment)
GA_E2_res <- hist(GA_E2_fragment$GA_E2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_E2_res$breaks),
     y = c(0, 0, GA_E2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_E2 Sample reads stat")

GA_E3_fragment <- GA_E3$count[GA_E3$count <= 5000]
head(GA_E2_fragment)
GA_E3_fragment <- data.frame(GA_E3_fragment)
GA_E3_res <- hist(GA_E3_fragment$GA_E3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_E3_res$breaks),
     y = c(0, 0, GA_E3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_E3 Sample reads stat")



GA <- cbind(GA_C1,GA_C2[,5],GA_C3[,5],GA_E1[,5],GA_E2[,5],GA_E3[,5])
GA_mtx <- as.data.frame(GA[,5:10])
row.names(GA_mtx) <- GA$id
GA_mtx <- na.omit(GA_mtx)
pheatmap(G6_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE) 
GAbaseMean = as.data.frame(colMeans(GA_mtx) ) 


G6_C1 <- fread("G6DD45A-C1.txt")
G6_C1[,"count"] <-G6_C1$meth +G6_C1$dmeth
G6_C2 <- fread("G6DD45A-C2.txt")
G6_C2[,"count"] <-G6_C2$meth +G6_C2$dmeth
G6_C3 <- fread("G6DD45A-C3.txt")
G6_C3[,"count"] <-G6_C3$meth +G6_C3$dmeth
G6_E1 <- fread("G6DD45A-OE1.txt")
G6_E1[,"count"] <-G6_E1$meth +G6_E1$dmeth
G6_E2 <- fread("G6DD45A-OE2.txt")
G6_E2[,"count"] <-G6_E2$meth +G6_E2$dmeth
G6_E3 <- fread("G6DD45A-OE3.txt")
G6_E3[,"count"] <-G6_E3$meth +G6_E3$dmeth

G6_C1_fragment <- G6_C1$count[G6_C1$count <= 5000]
head(G6_C1_fragment)
G6_C1_fragment <- data.frame(G6_C1_fragment)
G6_C1_res <- hist(G6_C1_fragment$G6_C1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_C1_res$breaks),
     y = c(0, 0, G6_C1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_C1 Sample reads stat")

G6_C2_fragment <- G6_C2$count[G6_C2$count <= 5000]
head(G6_C2_fragment)
G6_C2_fragment <- data.frame(G6_C2_fragment)
G6_C2_res <- hist(G6_C2_fragment$G6_C2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_C2_res$breaks),
     y = c(0, 0, G6_C2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_C2 Sample reads stat")

G6_C3_fragment <- G6_C3$count[G6_C3$count <= 5000]
head(G6_C2_fragment)
G6_C3_fragment <- data.frame(G6_C3_fragment)
G6_C3_res <- hist(G6_C3_fragment$G6_C3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_C3_res$breaks),
     y = c(0, 0, G6_C3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_C3 Sample reads stat")

G6_E1_fragment <- G6_E1$count[G6_E1$count <= 5000]
head(G6_E1_fragment)
G6_E1_fragment <- data.frame(G6_E1_fragment)
G6_E1_res <- hist(G6_E1_fragment$G6_E1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_E1_res$breaks),
     y = c(0, 0, G6_E1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_E1 Sample reads stat")

G6_E2_fragment <- G6_E2$count[G6_E2$count <= 5000]
head(G6_E2_fragment)
G6_E2_fragment <- data.frame(G6_E2_fragment)
G6_E2_res <- hist(G6_E2_fragment$G6_E2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_E2_res$breaks),
     y = c(0, 0, G6_E2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_E2 Sample reads stat")

G6_E3_fragment <- G6_E3$count[G6_E3$count <= 5000]
head(G6_E2_fragment)
G6_E3_fragment <- data.frame(G6_E3_fragment)
G6_E3_res <- hist(G6_E3_fragment$G6_E3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_E3_res$breaks),
     y = c(0, 0, G6_E3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_E3 Sample reads stat")


G6 <- cbind(G6_C1,G6_C2[,5],G6_C3[,5],G6_E1[,5],G6_E2[,5],G6_E3[,5])
G6_mtx <- as.data.frame(G6[,5:10])
row.names(G6_mtx) <- G6$id
G6_mtx <- na.omit(G6_mtx)
pheatmap(G6_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,
         angle_col = 45,treeheight_row = 0,
         border=FALSE) 
G6baseMean = as.data.frame(colMeans(G6_mtx) )


pmd_anno <- fread("global_pmd.bedGraph")
pmd_row_anno <- select(pmd_anno,V4,V9)
pmd_row_anno <- as.data.frame(pmd_row_anno[,-1])
row.names(pmd_row_anno) <- pmd_anno$V4
colnames(pmd_row_anno) <- "region_type"

pmd_regions <- row.names(subset(pmd_row_anno, subset= pmd_row_anno$region_type == "PMD"))
GAPMDMean = as.data.frame(colMeans(GA_mtx[pmd_regions,]) ) 
colnames(GAPMDMean) <- "GA-PMD"
colnames(GAbaseMean) <- "GA-base"



GA_C1_filt <- subset(GA_C1,subset = GA_C1$count >1000)
GA_C1_filt <- select(GA_C1_filt,id,'GADD45A-C1')
GA_C2_filt <- subset(GA_C2,subset = GA_C2$count >1000)
GA_C3_filt <- subset(GA_C3,subset = GA_C3$count >1000)
GA_E1_filt <- subset(GA_E1,subset = GA_E1$count >1000)
GA_E2_filt <- subset(GA_E2,subset = GA_E2$count >1000)
GA_E3_filt <- subset(GA_E3,subset = GA_E2$count >1000)
GA_C1_filt <- select(GA_C1_filt,id,'GADD45A-C1')
GA_C2_filt <- select(GA_C2_filt,id,'GADD45A-C2')
GA_C3_filt <- select(GA_C3_filt,id,'GADD45A-C3')
GA_E1_filt <- select(GA_E1_filt,id,'GADD45A-OE1')
GA_E2_filt <- select(GA_E2_filt,id,'GADD45A-OE2')
GA_E3_filt <- select(GA_E3_filt,id,'GADD45A-OE3')
GA_id <- intersect(GA_C1_filt$id,GA_C2_filt$id)
GA_id <- intersect(GA_id,GA_C3_filt$id)
GA_id <- intersect(GA_id,GA_E1_filt$id)
GA_id <- intersect(GA_id,GA_E2_filt$id)
GA_id <- intersect(GA_id,GA_E3_filt$id)
GA_C1_filt2 <- subset(GA_C1_filt, subset= GA_C1_filt$id %in% GA_id)
GA_C2_filt2 <- subset(GA_C2_filt, subset= GA_C2_filt$id %in% GA_id)
GA_C3_filt2 <- subset(GA_C3_filt, subset= GA_C3_filt$id %in% GA_id)
GA_E1_filt2 <- subset(GA_E1_filt, subset= GA_E1_filt$id %in% GA_id)
GA_E2_filt2 <- subset(GA_E2_filt, subset= GA_E2_filt$id %in% GA_id)
GA_E3_filt2 <- subset(GA_E3_filt, subset= GA_E3_filt$id %in% GA_id)
GA_filt <- cbind(GA_C1_filt2,GA_C2_filt2[,2],GA_C3_filt2[,2],GA_E1_filt2[,2],GA_E2_filt2[,2],GA_E3_filt2[,2])
GA_filt <- as.data.frame(GA_filt[,2:7])
row.names(GA_filt) <- GA_C1_filt2$id
GA_filt <- na.omit(GA_filt)
pheatmap(GA_filt,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,annotation_row = pmd_row_anno,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE) 
GAbaseMean = as.data.frame(colMeans(GA_filt) )



G6_C1_filt <- subset(G6_C1,subset = G6_C1$count >1000)
G6_C2_filt <- subset(G6_C2,subset = G6_C2$count >1000)
G6_C3_filt <- subset(G6_C3,subset = G6_C3$count >1000)
G6_E1_filt <- subset(G6_E1,subset = G6_E1$count >1000)
G6_E2_filt <- subset(G6_E2,subset = G6_E2$count >1000)
G6_E3_filt <- subset(G6_E3,subset = G6_E2$count >1000)
G6_C1_filt <- subset(G6_C1,subset = G6_C1$count >1000)
G6_C2_filt <- subset(G6_C2,subset = G6_C2$count >1000)
G6_C3_filt <- subset(G6_C3,subset = G6_C3$count >1000)
G6_E1_filt <- subset(G6_E1,subset = G6_E1$count >1000)
G6_E2_filt <- subset(G6_E2,subset = G6_E2$count >1000)
G6_E3_filt <- subset(G6_E3,subset = G6_E2$count >1000)
G6_C1_filt <- select(G6_C1_filt,id,'SNHG6-C1')
G6_C2_filt <- select(G6_C2_filt,id,'SNHG6-C2')
G6_C3_filt <- select(G6_C3_filt,id,'SNHG6-C3')
G6_E1_filt <- select(G6_E1_filt,id,'SNHG6-OE1')
G6_E2_filt <- select(G6_E2_filt,id,'SNHG6-OE2')
G6_E3_filt <- select(G6_E3_filt,id,'SNHG6-OE3')
G6_id <- intersect(G6_C1_filt$id,G6_C2_filt$id)
G6_id <- intersect(G6_id,G6_C3_filt$id)
G6_id <- intersect(G6_id,G6_E1_filt$id)
G6_id <- intersect(G6_id,G6_E2_filt$id)
G6_id <- intersect(G6_id,G6_E3_filt$id)
G6_C1_filt2 <- subset(G6_C1_filt, subset= G6_C1_filt$id %in% G6_id)
G6_C2_filt2 <- subset(G6_C2_filt, subset= G6_C2_filt$id %in% G6_id)
G6_C3_filt2 <- subset(G6_C3_filt, subset= G6_C3_filt$id %in% G6_id)
G6_E1_filt2 <- subset(G6_E1_filt, subset= G6_E1_filt$id %in% G6_id)
G6_E2_filt2 <- subset(G6_E2_filt, subset= G6_E2_filt$id %in% G6_id)
G6_E3_filt2 <- subset(G6_E3_filt, subset= G6_E3_filt$id %in% G6_id)
G6_filt <- cbind(G6_C1_filt2,G6_C2_filt2[,2],G6_C3_filt2[,2],G6_E1_filt2[,2],G6_E2_filt2[,2],G6_E3_filt2[,2])
G6_filt <- as.data.frame(G6_filt[,2:7])
row.names(G6_filt) <- G6_C1_filt2$id
G6_filt <- na.omit(G6_filt)
pheatmap(G6_filt,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,annotation_row = pmd_row_anno,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE) 
G6baseMean = as.data.frame(colMeans(G6_filt) )
G6PMDMean = as.data.frame(colMeans(G6_mtx[pmd_regions,]) ) 
colnames(G6PMDMean) <- "G6-PMD"
colnames(G6baseMean) <- "G6-base"


G6_test_base <- G6_C1_filt2 
G6_test_base[,"fc"] <- G6_test_base$pvalue
colnames(G6_test_base) <- c("id","pvalue","fc")
for (i in 1:nrow(G6_filt)) {
  test_kw<-data.frame(
    ratio=unlist(G6_filt[i,]),
    group=c(rep("control",3),
            rep("oe",3))
  )
  test_kw_result <- kruskal.test(ratio~group,data=test_kw)
  G6_test_base[i,2] <- test_kw_result$p.value
  G6_test_base[i,3] <- mean(unlist(test_kw[1:3,1]))/mean(unlist(test_kw[4:6,1]))
}

G6_test_fragment <- G6_test_base$pvalue

G6_test_fragment <- data.frame(G6_test_fragment)
G6_test_res <- hist(G6_test_fragment$G6_test_fragment, breaks = 5000, plot = FALSE)
plot(x = c(0, G6_test_res$breaks),
     y = c(0, 0, G6_test_res$counts),
     type = "l", col = "blue",
     xlab = "pvalue",
     ylab = "bin nums",
     main = "G6_test pvaluse stat")


GA_test_base <- GA_C1_filt2 
colnames(GA_test_base) <- c("id","pvalue")
GA_test_base[,"fc"] <- GA_test_base$pvalue
for (i in 1:nrow(GA_filt)) {
  test_kw<-data.frame(
    ratio=unlist(GA_filt[i,]),
    group=c(rep("control",3),
            rep("oe",3))
  )
  test_kw_result <- kruskal.test(ratio~group,data=test_kw)
  GA_test_base[i,2] <- test_kw_result$p.value
  GA_test_base[i,3] <- mean(unlist(test_kw[1:3,1]))/mean(unlist(test_kw[4:6,1]))
}



GA_test_sig <- subset(GA_test_base,subset = GA_test_base$pvalue < 0.05)
G6_test_sig <- subset(G6_test_base,subset = G6_test_base$pvalue < 0.05)
GA_test_sig_demeth <-  subset(GA_test_sig,subset = GA_test_sig$fc > 1)
GA_test_sig_admeth <-  subset(GA_test_sig,subset = GA_test_sig$fc < 1)
G6_test_sig_demeth <-  subset(G6_test_sig,subset = G6_test_sig$fc > 1)
G6_test_sig_admeth <-  subset(G6_test_sig,subset = G6_test_sig$fc < 1)

ga_dmr_regions <- as.data.frame(str_split_fixed(GA_test_sig_demeth$id,"_",2))
colnames(ga_dmr_regions) <- c("chr","start")
ga_dmr_regions$pos <- GA_test_sig_demeth$id
ga_dmr_regions$start <- as.numeric(ga_dmr_regions$start)
ga_dmr_regions$end <- ga_dmr_regions$start+1
ga_dmr_regions$start <- ga_dmr_regions$start*100000
ga_dmr_regions$end <- ga_dmr_regions$end*100000
ga_dmr_regions$end <- as.numeric(ga_dmr_regions$end)
ga_dmr_regions$start <- as.numeric(ga_dmr_regions$start)
ga_dmr_regions <- select(ga_dmr_regions,chr,start,end,pos)
write.table(ga_dmr_regions, "ga_dmr_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

ga_global_regions <- as.data.frame(str_split_fixed(GA_test_base$id,"_",2))
colnames(ga_global_regions) <- c("chr","start")
ga_global_regions$pos <- GA_test_base$id
ga_global_regions$start <- as.numeric(ga_global_regions$start)
ga_global_regions$end <- ga_global_regions$start+1
ga_global_regions$start <- ga_global_regions$start*100000
ga_global_regions$end <- ga_global_regions$end*100000
ga_global_regions$end <- as.numeric(ga_global_regions$end)
ga_global_regions$start <- as.numeric(ga_global_regions$start)
ga_global_regions <- select(ga_global_regions,chr,start,end,pos)
write.table(ga_global_regions, "ga_global_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

g6_dmr_regions <- as.data.frame(str_split_fixed(G6_test_sig_demeth$id,"_",2))
colnames(g6_dmr_regions) <- c("chr","start")
g6_dmr_regions$pos <- G6_test_sig_demeth$id
g6_dmr_regions$start <- as.numeric(g6_dmr_regions$start)
g6_dmr_regions$end <- g6_dmr_regions$start+1
g6_dmr_regions$start <- g6_dmr_regions$start*100000
g6_dmr_regions$end <- g6_dmr_regions$end*100000
g6_dmr_regions$end <- as.numeric(g6_dmr_regions$end)
g6_dmr_regions$start <- as.numeric(g6_dmr_regions$start)
g6_dmr_regions <- select(g6_dmr_regions,chr,start,end,pos)
write.table(g6_dmr_regions, "g6_dmr_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

g6_global_regions <- as.data.frame(str_split_fixed(G6_test_base$id,"_",2))
colnames(g6_global_regions) <- c("chr","start")
g6_global_regions$pos <- G6_test_base$id
g6_global_regions$start <- as.numeric(g6_global_regions$start)
g6_global_regions$end <- g6_global_regions$start+1
g6_global_regions$start <- g6_global_regions$start*100000
g6_global_regions$end <- g6_global_regions$end*100000
g6_global_regions$end <- as.numeric(g6_global_regions$end)
g6_global_regions$start <- as.numeric(g6_global_regions$start)
g6_global_regions <- select(g6_global_regions,chr,start,end,pos)
write.table(g6_global_regions, "g6_global_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

g6_admr_regions <- as.data.frame(str_split_fixed(G6_test_sig_admeth$id,"_",2))
colnames(g6_admr_regions) <- c("chr","start")
g6_admr_regions$pos <- G6_test_sig_admeth$id
g6_admr_regions$start <- as.numeric(g6_admr_regions$start)
g6_admr_regions$end <- g6_admr_regions$start+1
g6_admr_regions$start <- g6_admr_regions$start*100000
g6_admr_regions$end <- g6_admr_regions$end*100000
g6_admr_regions$end <- as.numeric(g6_admr_regions$end)
g6_admr_regions$start <- as.numeric(g6_admr_regions$start)
g6_admr_regions <- select(g6_admr_regions,chr,start,end,pos)
write.table(g6_admr_regions, "g6_admr_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


write.table(g6_admr_regions[,1:3], "g6_admr_regions2.bed",sep = "\t",quote=F,row.names = F,col.names = F)















GA_C1_bed <- as.data.frame(str_split_fixed(GA_C1_filt$id,"_",2))
colnames(GA_C1_bed) <- c("chr","start")
GA_C1_bed$cite <- GA_C1_filt$id
GA_C1_bed$meth_level <- GA_C1_filt$`GADD45A-C1`
GA_C1_bed$start <- as.numeric(GA_C1_bed$start)
GA_C1_bed$end <- GA_C1_bed$start+1
GA_C1_bed$start <- GA_C1_bed$start*100000
GA_C1_bed$end <- GA_C1_bed$end*100000
GA_C1_bed$end <- as.numeric(GA_C1_bed$end)
GA_C1_bed$start <- as.numeric(GA_C1_bed$start)
GA_C1_bed <- select(GA_C1_bed,chr,start,end,meth_level)
write.table(GA_C1_bed, "GA_C1.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



GA_C2_bed <- as.data.frame(str_split_fixed(GA_C2_filt$id,"_",2))
colnames(GA_C2_bed) <- c("chr","start")
GA_C2_bed$cite <- GA_C2_filt$id
GA_C2_bed$meth_level <- GA_C2_filt$`GADD45A-C2`
GA_C2_bed$start <- as.numeric(GA_C2_bed$start)
GA_C2_bed$end <- GA_C2_bed$start+1
GA_C2_bed$start <- GA_C2_bed$start*100000
GA_C2_bed$end <- GA_C2_bed$end*100000
GA_C2_bed$end <- as.numeric(GA_C2_bed$end)
GA_C2_bed$start <- as.numeric(GA_C2_bed$start)
GA_C2_bed <- select(GA_C2_bed,chr,start,end,meth_level)
write.table(GA_C2_bed, "GA_C2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)


GA_C3_bed <- as.data.frame(str_split_fixed(GA_C3_filt$id,"_",2))
colnames(GA_C3_bed) <- c("chr","start")
GA_C3_bed$cite <- GA_C3_filt$id
GA_C3_bed$meth_level <- GA_C3_filt$`GADD45A-C3`
GA_C3_bed$start <- as.numeric(GA_C3_bed$start)
GA_C3_bed$end <- GA_C3_bed$start+1
GA_C3_bed$start <- GA_C3_bed$start*100000
GA_C3_bed$end <- GA_C3_bed$end*100000
GA_C3_bed$end <- as.numeric(GA_C3_bed$end)
GA_C3_bed$start <- as.numeric(GA_C3_bed$start)
GA_C3_bed <- select(GA_C3_bed,chr,start,end,meth_level)
write.table(GA_C3_bed, "GA_C3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



GA_E1_bed <- as.data.frame(str_split_fixed(GA_E1_filt$id,"_",2))
colnames(GA_E1_bed) <- c("chr","start")
GA_E1_bed$cite <- GA_E1_filt$id
GA_E1_bed$meth_level <- GA_E1_filt$`GADD45A-OE1`
GA_E1_bed$start <- as.numeric(GA_E1_bed$start)
GA_E1_bed$end <- GA_E1_bed$start+1
GA_E1_bed$start <- GA_E1_bed$start*100000
GA_E1_bed$end <- GA_E1_bed$end*100000
GA_E1_bed$end <- as.numeric(GA_E1_bed$end)
GA_E1_bed$start <- as.numeric(GA_E1_bed$start)
GA_E1_bed <- select(GA_E1_bed,chr,start,end,meth_level)
write.table(GA_E1_bed, "GA_E1.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)




GA_E2_bed <- as.data.frame(str_split_fixed(GA_E2_filt$id,"_",2))
colnames(GA_E2_bed) <- c("chr","start")
GA_E2_bed$cite <- GA_E2_filt$id
GA_E2_bed$meth_level <- GA_E2_filt$`GADD45A-OE2`
GA_E2_bed$start <- as.numeric(GA_E2_bed$start)
GA_E2_bed$end <- GA_E2_bed$start+1
GA_E2_bed$start <- GA_E2_bed$start*100000
GA_E2_bed$end <- GA_E2_bed$end*100000
GA_E2_bed$end <- as.numeric(GA_E2_bed$end)
GA_E2_bed$start <- as.numeric(GA_E2_bed$start)
GA_E2_bed <- select(GA_E2_bed,chr,start,end,meth_level)
write.table(GA_E2_bed, "GA_E2.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)





GA_E3_bed <- as.data.frame(str_split_fixed(GA_E3_filt$id,"_",2))
colnames(GA_E3_bed) <- c("chr","start")
GA_E3_bed$cite <- GA_E3_filt$id
GA_E3_bed$meth_level <- GA_E3_filt$`GADD45A-OE3`
GA_E3_bed$start <- as.numeric(GA_E3_bed$start)
GA_E3_bed$end <- GA_E3_bed$start+1
GA_E3_bed$start <- GA_E3_bed$start*100000
GA_E3_bed$end <- GA_E3_bed$end*100000
GA_E3_bed$end <- as.numeric(GA_E3_bed$end)
GA_E3_bed$start <- as.numeric(GA_E3_bed$start)
GA_E3_bed <- select(GA_E3_bed,chr,start,end,meth_level)
write.table(GA_E3_bed, "GA_E3.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

