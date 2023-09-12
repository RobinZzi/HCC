library(VennDiagram)
rm(list=ls())
save.image("oe_10kbin.Rdata")
setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report/10kbin")

setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report/report")
fs <- list.files()
for (i in 1:length(fs)) {
  setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report/report")
  data <- fread(paste(fs[i]))
  data <- subset(data,data$V1 %in% chrs)
  colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
  data[,3] <- round(data[,2]/10000)
  data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
  data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
  sample_id <-str_split_fixed(fs[i],"_report.txt",2)[1]
  data[,"count"] <- data[,2]+data[,3]
  data[,paste(sample_id)] <- data[,2]/(data[,2]+data[,3])
  setwd("/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report/10kbin")
  write.table(data,paste(sample_id,"_10kbin.txt",sep = ""))
}




GA_C1 <- fread("GADD45A-C1_10kbin.txt")
GA_C2 <- fread("GADD45A-C2_10kbin.txt")
GA_C3 <- fread("GADD45A-C3_10kbin.txt")
GA_E1 <- fread("GADD45A-OE1_10kbin.txt")
GA_E2 <- fread("GADD45A-OE2_10kbin.txt")
GA_E3 <- fread("GADD45A-OE3_10kbin.txt")



GA_C1_fragment <- GA_C1$count[GA_C1$count <= 1000]
head(GA_C1_fragment)
GA_C1_fragment <- data.frame(GA_C1_fragment)
GA_C1_res <- hist(GA_C1_fragment$GA_C1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_C1_res$breaks),
     y = c(0, 0, GA_C1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_C1 Sample reads stat")

GA_C2_fragment <- GA_C2$count[GA_C2$count <= 1000]
head(GA_C2_fragment)
GA_C2_fragment <- data.frame(GA_C2_fragment)
GA_C2_res <- hist(GA_C2_fragment$GA_C2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_C2_res$breaks),
     y = c(0, 0, GA_C2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_C2 Sample reads stat")

GA_C3_fragment <- GA_C3$count[GA_C3$count <= 1000]
head(GA_C2_fragment)
GA_C3_fragment <- data.frame(GA_C3_fragment)
GA_C3_res <- hist(GA_C3_fragment$GA_C3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_C3_res$breaks),
     y = c(0, 0, GA_C3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_C3 Sample reads stat")

GA_E1_fragment <- GA_E1$count[GA_E1$count <= 1000]
head(GA_E1_fragment)
GA_E1_fragment <- data.frame(GA_E1_fragment)
GA_E1_res <- hist(GA_E1_fragment$GA_E1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_E1_res$breaks),
     y = c(0, 0, GA_E1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_E1 Sample reads stat")

GA_E2_fragment <- GA_E2$count[GA_E2$count <= 1000]
head(GA_E2_fragment)
GA_E2_fragment <- data.frame(GA_E2_fragment)
GA_E2_res <- hist(GA_E2_fragment$GA_E2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_E2_res$breaks),
     y = c(0, 0, GA_E2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_E2 Sample reads stat")

GA_E3_fragment <- GA_E3$count[GA_E3$count <= 1000]
head(GA_E2_fragment)
GA_E3_fragment <- data.frame(GA_E3_fragment)
GA_E3_res <- hist(GA_E3_fragment$GA_E3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_E3_res$breaks),
     y = c(0, 0, GA_E3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA_E3 Sample reads stat")



GA_C1_filt <- subset(GA_C1,subset = GA_C1$count >100)
GA_C2_filt <- subset(GA_C2,subset = GA_C2$count >100)
GA_C3_filt <- subset(GA_C3,subset = GA_C3$count >100)
GA_E1_filt <- subset(GA_E1,subset = GA_E1$count >100)
GA_E2_filt <- subset(GA_E2,subset = GA_E2$count >100)
GA_E3_filt <- subset(GA_E3,subset = GA_E2$count >100)
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
ga_pmd_regions <- ga_global_pmd$V4
GAPMDMean = as.data.frame(colMeans(GA_filt[ga_pmd_regions,]) ) 
colnames(GAPMDMean) <- "GA-PMD"
colnames(GAbaseMean) <- "GA-base"







if(T){
  GA_C1_bed <- as.data.frame(str_split_fixed(GA_C1_filt$id,"_",2))
  colnames(GA_C1_bed) <- c("chr","start")
  GA_C1_bed$cite <- GA_C1_filt$id
  GA_C1_bed$meth_level <- GA_C1_filt$`GADD45A-C1`
  GA_C1_bed$start <- as.numeric(GA_C1_bed$start)
  GA_C1_bed$end <- GA_C1_bed$start+1
  GA_C1_bed$start <- GA_C1_bed$start*10000
  GA_C1_bed$end <- GA_C1_bed$end*10000
  GA_C1_bed$end <- as.numeric(GA_C1_bed$end)
  GA_C1_bed$start <- as.numeric(GA_C1_bed$start)
  GA_C1_bed <- select(GA_C1_bed,chr,start,end,meth_level)
  write.table(GA_C1_bed, "GA_C1_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  GA_C2_bed <- as.data.frame(str_split_fixed(GA_C2_filt$id,"_",2))
  colnames(GA_C2_bed) <- c("chr","start")
  GA_C2_bed$cite <- GA_C2_filt$id
  GA_C2_bed$meth_level <- GA_C2_filt$`GADD45A-C2`
  GA_C2_bed$start <- as.numeric(GA_C2_bed$start)
  GA_C2_bed$end <- GA_C2_bed$start+1
  GA_C2_bed$start <- GA_C2_bed$start*10000
  GA_C2_bed$end <- GA_C2_bed$end*10000
  GA_C2_bed$end <- as.numeric(GA_C2_bed$end)
  GA_C2_bed$start <- as.numeric(GA_C2_bed$start)
  GA_C2_bed <- select(GA_C2_bed,chr,start,end,meth_level)
  write.table(GA_C2_bed, "GA_C2_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  GA_C3_bed <- as.data.frame(str_split_fixed(GA_C3_filt$id,"_",2))
  colnames(GA_C3_bed) <- c("chr","start")
  GA_C3_bed$cite <- GA_C3_filt$id
  GA_C3_bed$meth_level <- GA_C3_filt$`GADD45A-C3`
  GA_C3_bed$start <- as.numeric(GA_C3_bed$start)
  GA_C3_bed$end <- GA_C3_bed$start+1
  GA_C3_bed$start <- GA_C3_bed$start*10000
  GA_C3_bed$end <- GA_C3_bed$end*10000
  GA_C3_bed$end <- as.numeric(GA_C3_bed$end)
  GA_C3_bed$start <- as.numeric(GA_C3_bed$start)
  GA_C3_bed <- select(GA_C3_bed,chr,start,end,meth_level)
  write.table(GA_C3_bed, "GA_C3_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  GA_E1_bed <- as.data.frame(str_split_fixed(GA_E1_filt$id,"_",2))
  colnames(GA_E1_bed) <- c("chr","start")
  GA_E1_bed$cite <- GA_E1_filt$id
  GA_E1_bed$meth_level <- GA_E1_filt$`GADD45A-OE1`
  GA_E1_bed$start <- as.numeric(GA_E1_bed$start)
  GA_E1_bed$end <- GA_E1_bed$start+1
  GA_E1_bed$start <- GA_E1_bed$start*10000
  GA_E1_bed$end <- GA_E1_bed$end*10000
  GA_E1_bed$end <- as.numeric(GA_E1_bed$end)
  GA_E1_bed$start <- as.numeric(GA_E1_bed$start)
  GA_E1_bed <- select(GA_E1_bed,chr,start,end,meth_level)
  write.table(GA_E1_bed, "GA_E1_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  
  GA_E2_bed <- as.data.frame(str_split_fixed(GA_E2_filt$id,"_",2))
  colnames(GA_E2_bed) <- c("chr","start")
  GA_E2_bed$cite <- GA_E2_filt$id
  GA_E2_bed$meth_level <- GA_E2_filt$`GADD45A-OE2`
  GA_E2_bed$start <- as.numeric(GA_E2_bed$start)
  GA_E2_bed$end <- GA_E2_bed$start+1
  GA_E2_bed$start <- GA_E2_bed$start*10000
  GA_E2_bed$end <- GA_E2_bed$end*10000
  GA_E2_bed$end <- as.numeric(GA_E2_bed$end)
  GA_E2_bed$start <- as.numeric(GA_E2_bed$start)
  GA_E2_bed <- select(GA_E2_bed,chr,start,end,meth_level)
  write.table(GA_E2_bed, "GA_E2_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  
  
  GA_E3_bed <- as.data.frame(str_split_fixed(GA_E3_filt$id,"_",2))
  colnames(GA_E3_bed) <- c("chr","start")
  GA_E3_bed$cite <- GA_E3_filt$id
  GA_E3_bed$meth_level <- GA_E3_filt$`GADD45A-OE3`
  GA_E3_bed$start <- as.numeric(GA_E3_bed$start)
  GA_E3_bed$end <- GA_E3_bed$start+1
  GA_E3_bed$start <- GA_E3_bed$start*10000
  GA_E3_bed$end <- GA_E3_bed$end*10000
  GA_E3_bed$end <- as.numeric(GA_E3_bed$end)
  GA_E3_bed$start <- as.numeric(GA_E3_bed$start)
  GA_E3_bed <- select(GA_E3_bed,chr,start,end,meth_level)
  write.table(GA_E3_bed, "GA_E3_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
}



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
  GA_test_base[i,3] <- mean(unlist(test_kw[1:3,1]))-mean(unlist(test_kw[4:6,1]))
}



GA_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[1:3], row[4:6])$p.value
})
GA_dt <- apply(GA_filt, 1, function(row) {
  mean(row[1:3])-mean(row[4:6])
})

GA_adjusted_p_values <- p.adjust(GA_p_values, method  = "fdr")
GA_adjusted_p_values_b <- p.adjust(GA_p_values, method  = "bonferroni")
GA_test <- cbind(GA_p_values,GA_adjusted_p_values,GA_adjusted_p_values_b,GA_dt)
GA_test <- as.data.frame(GA_test)
colnames(GA_test) <- c("pvalue","bf_adj","fdr","dt")
GA_test_significant <- subset(GA_test,subset =  pvalue< 0.05)
GA_ttest_sig_demeth <- subset(GA_test_significant , subset = dt>0)
GA_ttest_sig_admeth <- subset(GA_test_significant , subset = dt<0)


GA_test_sig <- subset(GA_test_base,subset = GA_test_base$pvalue < 0.05)
GA_test_sig_demeth <-  subset(GA_test_sig,subset = GA_test_sig$fc > 0)
GA_test_sig_admeth <-  subset(GA_test_sig,subset = GA_test_sig$fc < 0)

ga_demeth_Venn <- list(ttest = row.names(GA_ttest_sig_demeth), kwtest = GA_test_sig_demeth$id)
ga_admeth_Venn <- list(ttest = row.names(GA_ttest_sig_admeth), kwtest = GA_test_sig_admeth$id)


venn.diagram(ga_demeth_Venn, filename = 'ga_demeth_Venn.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB'), cex = 1.5, fontfamily = 'serif')

venn.diagram(ga_admeth_Venn, filename = 'ga_admeth_Venn.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB'), cex = 1.5, fontfamily = 'serif')


wtest_demeth <- intersect(GA_test_sig_demeth$id,row.names(GA_ttest_sig_demeth))
wtest_admeth <- intersect(GA_test_sig_admeth$id,row.names(GA_ttest_sig_admeth))

if(T){
ga_dmr_regions <- as.data.frame(str_split_fixed(wtest_demeth,"_",2))
colnames(ga_dmr_regions) <- c("chr","start")
ga_dmr_regions$pos <- wtest_demeth
ga_dmr_regions$start <- as.numeric(ga_dmr_regions$start)
ga_dmr_regions$end <- ga_dmr_regions$start+1
ga_dmr_regions$start <- ga_dmr_regions$start*10000
ga_dmr_regions$end <- ga_dmr_regions$end*10000
ga_dmr_regions$end <- as.numeric(ga_dmr_regions$end)
ga_dmr_regions$start <- as.numeric(ga_dmr_regions$start)
ga_dmr_regions <- select(ga_dmr_regions,chr,start,end,pos)
write.table(ga_dmr_regions, "ga_10k_dmr_regions.bed",sep = "\t",quote=F,row.names = F,col.names = F)



ga_admr_regions <- as.data.frame(str_split_fixed(wtest_admeth,"_",2))
colnames(ga_admr_regions) <- c("chr","start")
ga_admr_regions$pos <- wtest_admeth
ga_admr_regions$start <- as.numeric(ga_admr_regions$start)
ga_admr_regions$end <- ga_admr_regions$start+1
ga_admr_regions$start <- ga_admr_regions$start*10000
ga_admr_regions$end <- ga_admr_regions$end*10000
ga_admr_regions$end <- as.numeric(ga_admr_regions$end)
ga_admr_regions$start <- as.numeric(ga_admr_regions$start)
ga_admr_regions <- select(ga_admr_regions,chr,start,end,pos)
write.table(ga_admr_regions, "ga_10k_admr_regions.bed",sep = "\t",quote=F,row.names = F,col.names = F)






ga_global_regions <- as.data.frame(str_split_fixed(row.names(GA_filt),"_",2))
colnames(ga_global_regions ) <- c("chr","start")
ga_global_regions$pos <- row.names(GA_filt)
ga_global_regions$start <- as.numeric(ga_global_regions$start)
ga_global_regions$end <- ga_global_regions$start+1
ga_global_regions$start <- ga_global_regions$start*10000
ga_global_regions$end <- ga_global_regions$end*10000
ga_global_regions$end <- as.numeric(ga_global_regions$end)
ga_global_regions$start <- as.numeric(ga_global_regions$start)
ga_global_regions <- select(ga_global_regions,chr,start,end,pos)
write.table(ga_global_regions, "ga_10k_global_regions.bed",sep = "\t",quote=F,row.names = F,col.names = F)
}





if(T){
  
  
  ga_global_H3K9me3 <- fread("ga_global_H3K9me3.bedGraph")
  ga_global_H3K9me3_count <- select(ga_global_H3K9me3,V14,V15)
  ga_global_H3K9me3_count <- aggregate(ga_global_H3K9me3_count[,2],by=list(id=ga_global_H3K9me3_count$V14),FUN=sum)
  ga_global_H3K9me3_count$hm <- 'H3K9me3'
  ga_global_H3K9me3_count$rg <- 'ga_global'
  ga_global_H3K9me3_count$V15 <- ga_global_H3K9me3_count$V15/nrow(ga_global_regions)
  
  
  ga_dmr_H3K9me3 <- fread("ga_dmr_H3K9me3.bedGraph")
  ga_dmr_H3K9me3_count <- select(ga_dmr_H3K9me3,V14,V15)
  ga_dmr_H3K9me3_count <- aggregate(ga_dmr_H3K9me3_count[,2],by=list(id=ga_dmr_H3K9me3_count$V14),FUN=sum)
  ga_dmr_H3K9me3_count$hm <- 'H3K9me3'
  ga_dmr_H3K9me3_count$rg <- 'ga_dmr'
  ga_dmr_H3K9me3_count$V15 <- ga_dmr_H3K9me3_count$V15/nrow(ga_dmr_regions)
  
  
  ga_admr_H3K9me3 <- fread("ga_admr_H3K9me3.bedGraph")
  ga_admr_H3K9me3_count <- select(ga_admr_H3K9me3,V14,V15)
  ga_admr_H3K9me3_count <- aggregate(ga_admr_H3K9me3_count[,2],by=list(id=ga_admr_H3K9me3_count$V14),FUN=sum)
  ga_admr_H3K9me3_count$hm <- 'H3K9me3'
  ga_admr_H3K9me3_count$rg <- 'ga_admr'
  ga_admr_H3K9me3_count$V15 <- ga_admr_H3K9me3_count$V15/nrow(ga_admr_regions)
  ###########################
  ga_global_H3K36me3 <- fread("ga_global_H3K36me3.bedGraph")
  ga_global_H3K36me3_count <- select(ga_global_H3K36me3,V14,V15)
  ga_global_H3K36me3_count <- aggregate(ga_global_H3K36me3_count[,2],by=list(id=ga_global_H3K36me3_count$V14),FUN=sum)
  ga_global_H3K36me3_count$hm <- 'H3K36me3'
  ga_global_H3K36me3_count$rg <- 'ga_global'
  ga_global_H3K36me3_count$V15 <- ga_global_H3K36me3_count$V15/nrow(ga_global_regions)
  
  ga_dmr_H3K36me3 <- fread("ga_dmr_H3K36me3.bedGraph")
  ga_dmr_H3K36me3_count <- select(ga_dmr_H3K36me3,V14,V15)
  ga_dmr_H3K36me3_count <- aggregate(ga_dmr_H3K36me3_count[,2],by=list(id=ga_dmr_H3K36me3_count$V14),FUN=sum)
  ga_dmr_H3K36me3_count$hm <- 'H3K36me3'
  ga_dmr_H3K36me3_count$rg <- 'ga_dmr'
  ga_dmr_H3K36me3_count$V15 <- ga_dmr_H3K36me3_count$V15/nrow(ga_dmr_regions)
  
  
  ga_admr_H3K36me3 <- fread("ga_admr_H3K36me3.bedGraph")
  ga_admr_H3K36me3_count <- select(ga_admr_H3K36me3,V14,V15)
  ga_admr_H3K36me3_count <- aggregate(ga_admr_H3K36me3_count[,2],by=list(id=ga_admr_H3K36me3_count$V14),FUN=sum)
  ga_admr_H3K36me3_count$hm <- 'H3K36me3'
  ga_admr_H3K36me3_count$rg <- 'ga_admr'
  ga_admr_H3K36me3_count$V15 <- ga_admr_H3K36me3_count$V15/nrow(ga_admr_regions)
  
  ###########################
  ga_global_H3K27me3 <- fread("ga_global_H3K27me3.bedGraph")
  ga_global_H3K27me3_count <- select(ga_global_H3K27me3,V14,V15)
  ga_global_H3K27me3_count <- aggregate(ga_global_H3K27me3_count[,2],by=list(id=ga_global_H3K27me3_count$V14),FUN=sum)
  ga_global_H3K27me3_count$hm <- 'H3K27me3'
  ga_global_H3K27me3_count$rg <- 'ga_global'
  ga_global_H3K27me3_count$V15 <- ga_global_H3K27me3_count$V15/nrow(ga_global_regions)
  
  
  ga_dmr_H3K27me3 <- fread("ga_dmr_H3K27me3.bedGraph")
  ga_dmr_H3K27me3_count <- select(ga_dmr_H3K27me3,V14,V15)
  ga_dmr_H3K27me3_count <- aggregate(ga_dmr_H3K27me3_count[,2],by=list(id=ga_dmr_H3K27me3_count$V14),FUN=sum)
  ga_dmr_H3K27me3_count$hm <- 'H3K27me3'
  ga_dmr_H3K27me3_count$rg <- 'ga_dmr'
  ga_dmr_H3K27me3_count$V15 <- ga_dmr_H3K27me3_count$V15/nrow(ga_dmr_regions)
  
  ga_admr_H3K27me3 <- fread("ga_admr_H3K27me3.bedGraph")
  ga_admr_H3K27me3_count <- select(ga_admr_H3K27me3,V14,V15)
  ga_admr_H3K27me3_count <- aggregate(ga_admr_H3K27me3_count[,2],by=list(id=ga_admr_H3K27me3_count$V14),FUN=sum)
  ga_admr_H3K27me3_count$hm <- 'H3K27me3'
  ga_admr_H3K27me3_count$rg <- 'ga_admr'
  ga_admr_H3K27me3_count$V15 <- ga_admr_H3K27me3_count$V15/nrow(ga_admr_regions)
  
  ###########################
  ga_global_H3K27ac <- fread("ga_global_H3K27ac.bedGraph")
  ga_global_H3K27ac_count <- select(ga_global_H3K27ac,V14,V15)
  ga_global_H3K27ac_count <- aggregate(ga_global_H3K27ac_count[,2],by=list(id=ga_global_H3K27ac_count$V14),FUN=sum)
  ga_global_H3K27ac_count$hm <- 'H3K27ac'
  ga_global_H3K27ac_count$rg <- 'ga_global'
  ga_global_H3K27ac_count$V15 <- ga_global_H3K27ac_count$V15/nrow(ga_global_regions)
  
  
  ga_dmr_H3K27ac <- fread("ga_dmr_H3K27ac.bedGraph")
  ga_dmr_H3K27ac_count <- select(ga_dmr_H3K27ac,V14,V15)
  ga_dmr_H3K27ac_count <- aggregate(ga_dmr_H3K27ac_count[,2],by=list(id=ga_dmr_H3K27ac_count$V14),FUN=sum)
  ga_dmr_H3K27ac_count$hm <- 'H3K27ac'
  ga_dmr_H3K27ac_count$rg <- 'ga_dmr'
  ga_dmr_H3K27ac_count$V15 <- ga_dmr_H3K27ac_count$V15/nrow(ga_dmr_regions)
  
  
  ga_admr_H3K27ac <- fread("ga_admr_H3K27ac.bedGraph")
  ga_admr_H3K27ac_count <- select(ga_admr_H3K27ac,V14,V15)
  ga_admr_H3K27ac_count <- aggregate(ga_admr_H3K27ac_count[,2],by=list(id=ga_admr_H3K27ac_count$V14),FUN=sum)
  ga_admr_H3K27ac_count$hm <- 'H3K27ac'
  ga_admr_H3K27ac_count$rg <- 'ga_admr'
  ga_admr_H3K27ac_count$V15 <- ga_admr_H3K27ac_count$V15/nrow(ga_admr_regions)
}
ga_hm_combine <- rbind(ga_global_H3K9me3_count,ga_dmr_H3K9me3_count,ga_global_H3K36me3_count,ga_dmr_H3K36me3_count,
                       ga_global_H3K27me3_count,ga_dmr_H3K27me3_count,
                       ga_global_H3K27ac_count,ga_dmr_H3K27ac_count,ga_admr_H3K27ac_count,ga_admr_H3K36me3_count,ga_admr_H3K9me3_count,ga_admr_H3K27me3_count)
colnames(ga_hm_combine) <- c("region","length","hm","region_type")
ggplot(ga_hm_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))



ggplot(ga_hm_combine)+
  geom_bar(aes(x=region_type,y=length,fill=hm),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))




ga_global_anno <- fread("ga_global_PMD.bedGraph")
ga_global_pmd <- subset(ga_global_anno,ga_global_anno$V9 == 'PMD')
ga_global_hmd <- subset(ga_global_anno,ga_global_anno$V9 == 'HMD')



ga_global_pmd_count <- select(ga_global_pmd,V4,V9,V11)
ga_global_pmd_count <- aggregate(ga_global_pmd_count[,3],by=list(id=ga_global_pmd_count$V9),FUN=sum)
ga_global_pmd_count$rg <- 'ga_global'
colnames(ga_global_pmd_count) <- c("hm","length","rg")
ga_global_pmd_count$length <- ga_global_pmd_count$length/nrow(ga_global_regions)

ga_global_hmd_count <- select(ga_global_hmd,V4,V9,V11)
ga_global_hmd_count <- aggregate(ga_global_hmd_count[,3],by=list(id=ga_global_hmd_count$V9),FUN=sum)
ga_global_hmd_count$rg <- 'ga_global'
colnames(ga_global_hmd_count) <-c("hm","length","rg")
ga_global_hmd_count$length <- ga_global_hmd_count$length/nrow(ga_global_regions)






ga_dmr_anno <- fread("ga_dmr_PMD.bedGraph")
ga_dmr_pmd <- subset(ga_dmr_anno,ga_dmr_anno$V9 == 'PMD')
ga_dmr_hmd <- subset(ga_dmr_anno,ga_dmr_anno$V9 == 'HMD')



ga_dmr_pmd_count <- select(ga_dmr_pmd,V4,V9,V11)
ga_dmr_pmd_count <- aggregate(ga_dmr_pmd_count[,3],by=list(id=ga_dmr_pmd_count$V9),FUN=sum)
ga_dmr_pmd_count$rg <- 'ga_dmr'
colnames(ga_dmr_pmd_count) <- c("hm","length","rg")
ga_dmr_pmd_count$length <- ga_dmr_pmd_count$length/nrow(ga_dmr_regions)

ga_dmr_hmd_count <- select(ga_dmr_hmd,V4,V9,V11)
ga_dmr_hmd_count <- aggregate(ga_dmr_hmd_count[,3],by=list(id=ga_dmr_hmd_count$V9),FUN=sum)
ga_dmr_hmd_count$rg <- 'ga_dmr'
colnames(ga_dmr_hmd_count) <- c("hm","length","rg")
ga_dmr_hmd_count$length <- ga_dmr_hmd_count$length/nrow(ga_dmr_regions)







ga_admr_anno <- fread("ga_admr_PMD.bedGraph")
ga_admr_pmd <- subset(ga_admr_anno,ga_admr_anno$V9 == 'PMD')
ga_admr_hmd <- subset(ga_admr_anno,ga_admr_anno$V9 == 'HMD')



ga_admr_pmd_count <- select(ga_admr_pmd,V4,V9,V11)
ga_admr_pmd_count <- aggregate(ga_admr_pmd_count[,3],by=list(id=ga_admr_pmd_count$V9),FUN=sum)
ga_admr_pmd_count$rg <- 'ga_admr'
colnames(ga_admr_pmd_count) <- c("hm","length","rg")
ga_admr_pmd_count$length <- ga_admr_pmd_count$length/nrow(ga_admr_regions)

ga_admr_hmd_count <- select(ga_admr_hmd,V4,V9,V11)
ga_admr_hmd_count <- aggregate(ga_admr_hmd_count[,3],by=list(id=ga_admr_hmd_count$V9),FUN=sum)
ga_admr_hmd_count$rg <- 'ga_admr'
colnames(ga_admr_hmd_count) <- c("hm","length","rg")
ga_admr_hmd_count$length <- ga_admr_hmd_count$length/nrow(ga_admr_regions)




ga_md_combine <- rbind(ga_admr_hmd_count,ga_admr_pmd_count,ga_dmr_hmd_count,ga_dmr_pmd_count,ga_global_hmd_count,ga_global_pmd_count)
colnames(ga_md_combine) <- c("hm","length","region_type")
ggplot(ga_md_combine)+
  geom_bar(aes(x=region_type,y=length,fill=hm),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))
ggplot(ga_md_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))





G6_C1 <- fread("SNHG6-C1_10kbin.txt")
G6_C2 <- fread("SNHG6-C2_10kbin.txt")
G6_C3 <- fread("SNHG6-C3_10kbin.txt")
G6_E1 <- fread("SNHG6-OE1_10kbin.txt")
G6_E2 <- fread("SNHG6-OE2_10kbin.txt")
G6_E3 <- fread("SNHG6-OE3_10kbin.txt")

G6_C1_fragment <- G6_C1$count[G6_C1$count <= 1000]
head(G6_C1_fragment)
G6_C1_fragment <- data.frame(G6_C1_fragment)
G6_C1_res <- hist(G6_C1_fragment$G6_C1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_C1_res$breaks),
     y = c(0, 0, G6_C1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_C1 Sample reads stat")

G6_C2_fragment <- G6_C2$count[G6_C2$count <= 1000]
head(G6_C2_fragment)
G6_C2_fragment <- data.frame(G6_C2_fragment)
G6_C2_res <- hist(G6_C2_fragment$G6_C2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_C2_res$breaks),
     y = c(0, 0, G6_C2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_C2 Sample reads stat")

G6_C3_fragment <- G6_C3$count[G6_C3$count <= 1000]
head(G6_C2_fragment)
G6_C3_fragment <- data.frame(G6_C3_fragment)
G6_C3_res <- hist(G6_C3_fragment$G6_C3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_C3_res$breaks),
     y = c(0, 0, G6_C3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_C3 Sample reads stat")

G6_E1_fragment <- G6_E1$count[G6_E1$count <= 1000]
head(G6_E1_fragment)
G6_E1_fragment <- data.frame(G6_E1_fragment)
G6_E1_res <- hist(G6_E1_fragment$G6_E1_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_E1_res$breaks),
     y = c(0, 0, G6_E1_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_E1 Sample reads stat")

G6_E2_fragment <- G6_E2$count[G6_E2$count <= 1000]
head(G6_E2_fragment)
G6_E2_fragment <- data.frame(G6_E2_fragment)
G6_E2_res <- hist(G6_E2_fragment$G6_E2_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_E2_res$breaks),
     y = c(0, 0, G6_E2_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_E2 Sample reads stat")

G6_E3_fragment <- G6_E3$count[G6_E3$count <= 1000]
head(G6_E2_fragment)
G6_E3_fragment <- data.frame(G6_E3_fragment)
G6_E3_res <- hist(G6_E3_fragment$G6_E3_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_E3_res$breaks),
     y = c(0, 0, G6_E3_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6_E3 Sample reads stat")



G6_C1_filt <- subset(G6_C1,subset = G6_C1$count >100)
G6_C2_filt <- subset(G6_C2,subset = G6_C2$count >100)
G6_C3_filt <- subset(G6_C3,subset = G6_C3$count >100)
G6_E1_filt <- subset(G6_E1,subset = G6_E1$count >100)
G6_E2_filt <- subset(G6_E2,subset = G6_E2$count >100)
G6_E3_filt <- subset(G6_E3,subset = G6_E2$count >100)
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








if(T){
  G6_C1_bed <- as.data.frame(str_split_fixed(G6_C1_filt$id,"_",2))
  colnames(G6_C1_bed) <- c("chr","start")
  G6_C1_bed$cite <- G6_C1_filt$id
  G6_C1_bed$meth_level <- G6_C1_filt$`SNHG6-C1`
  G6_C1_bed$start <- as.numeric(G6_C1_bed$start)
  G6_C1_bed$end <- G6_C1_bed$start+1
  G6_C1_bed$start <- G6_C1_bed$start*10000
  G6_C1_bed$end <- G6_C1_bed$end*10000
  G6_C1_bed$end <- as.numeric(G6_C1_bed$end)
  G6_C1_bed$start <- as.numeric(G6_C1_bed$start)
  G6_C1_bed <- select(G6_C1_bed,chr,start,end,meth_level)
  write.table(G6_C1_bed, "G6_C1_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  G6_C2_bed <- as.data.frame(str_split_fixed(G6_C2_filt$id,"_",2))
  colnames(G6_C2_bed) <- c("chr","start")
  G6_C2_bed$cite <- G6_C2_filt$id
  G6_C2_bed$meth_level <- G6_C2_filt$`SNHG6-C2`
  G6_C2_bed$start <- as.numeric(G6_C2_bed$start)
  G6_C2_bed$end <- G6_C2_bed$start+1
  G6_C2_bed$start <- G6_C2_bed$start*10000
  G6_C2_bed$end <- G6_C2_bed$end*10000
  G6_C2_bed$end <- as.numeric(G6_C2_bed$end)
  G6_C2_bed$start <- as.numeric(G6_C2_bed$start)
  G6_C2_bed <- select(G6_C2_bed,chr,start,end,meth_level)
  write.table(G6_C2_bed, "G6_C2_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  G6_C3_bed <- as.data.frame(str_split_fixed(G6_C3_filt$id,"_",2))
  colnames(G6_C3_bed) <- c("chr","start")
  G6_C3_bed$cite <- G6_C3_filt$id
  G6_C3_bed$meth_level <- G6_C3_filt$`SNHG6-C3`
  G6_C3_bed$start <- as.numeric(G6_C3_bed$start)
  G6_C3_bed$end <- G6_C3_bed$start+1
  G6_C3_bed$start <- G6_C3_bed$start*10000
  G6_C3_bed$end <- G6_C3_bed$end*10000
  G6_C3_bed$end <- as.numeric(G6_C3_bed$end)
  G6_C3_bed$start <- as.numeric(G6_C3_bed$start)
  G6_C3_bed <- select(G6_C3_bed,chr,start,end,meth_level)
  write.table(G6_C3_bed, "G6_C3_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  G6_E1_bed <- as.data.frame(str_split_fixed(G6_E1_filt$id,"_",2))
  colnames(G6_E1_bed) <- c("chr","start")
  G6_E1_bed$cite <- G6_E1_filt$id
  G6_E1_bed$meth_level <- G6_E1_filt$`SNHG6-OE1`
  G6_E1_bed$start <- as.numeric(G6_E1_bed$start)
  G6_E1_bed$end <- G6_E1_bed$start+1
  G6_E1_bed$start <- G6_E1_bed$start*10000
  G6_E1_bed$end <- G6_E1_bed$end*10000
  G6_E1_bed$end <- as.numeric(G6_E1_bed$end)
  G6_E1_bed$start <- as.numeric(G6_E1_bed$start)
  G6_E1_bed <- select(G6_E1_bed,chr,start,end,meth_level)
  write.table(G6_E1_bed, "G6_E1_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  
  G6_E2_bed <- as.data.frame(str_split_fixed(G6_E2_filt$id,"_",2))
  colnames(G6_E2_bed) <- c("chr","start")
  G6_E2_bed$cite <- G6_E2_filt$id
  G6_E2_bed$meth_level <- G6_E2_filt$`SNHG6-OE2`
  G6_E2_bed$start <- as.numeric(G6_E2_bed$start)
  G6_E2_bed$end <- G6_E2_bed$start+1
  G6_E2_bed$start <- G6_E2_bed$start*10000
  G6_E2_bed$end <- G6_E2_bed$end*10000
  G6_E2_bed$end <- as.numeric(G6_E2_bed$end)
  G6_E2_bed$start <- as.numeric(G6_E2_bed$start)
  G6_E2_bed <- select(G6_E2_bed,chr,start,end,meth_level)
  write.table(G6_E2_bed, "G6_E2_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  
  
  G6_E3_bed <- as.data.frame(str_split_fixed(G6_E3_filt$id,"_",2))
  colnames(G6_E3_bed) <- c("chr","start")
  G6_E3_bed$cite <- G6_E3_filt$id
  G6_E3_bed$meth_level <- G6_E3_filt$`SNHG6-OE3`
  G6_E3_bed$start <- as.numeric(G6_E3_bed$start)
  G6_E3_bed$end <- G6_E3_bed$start+1
  G6_E3_bed$start <- G6_E3_bed$start*10000
  G6_E3_bed$end <- G6_E3_bed$end*10000
  G6_E3_bed$end <- as.numeric(G6_E3_bed$end)
  G6_E3_bed$start <- as.numeric(G6_E3_bed$start)
  G6_E3_bed <- select(G6_E3_bed,chr,start,end,meth_level)
  G6_E3_bed <- na.omit(G6_E3_bed)
  write.table(G6_E3_bed, "G6_E3_10k.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
}



G6_test_base <- G6_C1_filt2 
colnames(G6_test_base) <- c("id","pvalue")
G6_test_base[,"fc"] <- G6_test_base$pvalue
for (i in 1:nrow(G6_filt)) {
  test_kw<-data.frame(
    ratio=unlist(G6_filt[i,]),
    group=c(rep("control",3),
            rep("oe",3))
  )
  test_kw_result <- kruskal.test(ratio~group,data=test_kw)
  G6_test_base[i,2] <- test_kw_result$p.value
  G6_test_base[i,3] <- mean(unlist(test_kw[1:3,1]))-mean(unlist(test_kw[4:6,1]))
}



G6_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[1:3], row[4:6])$p.value
})
G6_dt <- apply(G6_filt, 1, function(row) {
  mean(row[1:3])-mean(row[4:6])
})

G6_adjusted_p_values <- p.adjust(G6_p_values, method  = "fdr")
G6_adjusted_p_values_b <- p.adjust(G6_p_values, method  = "bonferroni")
G6_test <- cbind(G6_p_values,G6_adjusted_p_values,G6_adjusted_p_values_b,G6_dt)
G6_test <- as.data.frame(G6_test)
colnames(G6_test) <- c("pvalue","bf_adj","fdr","dt")
G6_test_significant <- subset(G6_test,subset =  pvalue< 0.05)
G6_ttest_sig_demeth <- subset(G6_test_significant , subset = dt>0)
G6_ttest_sig_admeth <- subset(G6_test_significant , subset = dt<0)


G6_test_sig <- subset(G6_test_base,subset = G6_test_base$pvalue < 0.05)
G6_test_sig_demeth <-  subset(G6_test_sig,subset = G6_test_sig$fc > 0)
G6_test_sig_admeth <-  subset(G6_test_sig,subset = G6_test_sig$fc < 0)

G6_demeth_Venn <- list(ttest = row.names(G6_ttest_sig_demeth), kwtest = G6_test_sig_demeth$id)
G6_admeth_Venn <- list(ttest = row.names(G6_ttest_sig_admeth), kwtest = G6_test_sig_admeth$id)


venn.diagram(G6_demeth_Venn, filename = 'G6_demeth_Venn.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB'), cex = 1.5, fontfamily = 'serif')

venn.diagram(G6_admeth_Venn, filename = 'G6_admeth_Venn.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB'), cex = 1.5, fontfamily = 'serif')


wtest_demeth <- intersect(G6_test_sig_demeth$id,row.names(G6_ttest_sig_demeth))
wtest_admeth <- intersect(G6_test_sig_admeth$id,row.names(G6_ttest_sig_admeth))

if(T){
  G6_dmr_regions <- as.data.frame(str_split_fixed(wtest_demeth,"_",2))
  colnames(G6_dmr_regions) <- c("chr","start")
  G6_dmr_regions$pos <- wtest_demeth
  G6_dmr_regions$start <- as.numeric(G6_dmr_regions$start)
  G6_dmr_regions$end <- G6_dmr_regions$start+1
  G6_dmr_regions$start <- G6_dmr_regions$start*10000
  G6_dmr_regions$end <- G6_dmr_regions$end*10000
  G6_dmr_regions$end <- as.numeric(G6_dmr_regions$end)
  G6_dmr_regions$start <- as.numeric(G6_dmr_regions$start)
  G6_dmr_regions <- select(G6_dmr_regions,chr,start,end,pos)
  write.table(G6_dmr_regions, "G6_10k_dmr_regions.bed",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  G6_admr_regions <- as.data.frame(str_split_fixed(wtest_admeth,"_",2))
  colnames(G6_admr_regions) <- c("chr","start")
  G6_admr_regions$pos <- wtest_admeth
  G6_admr_regions$start <- as.numeric(G6_admr_regions$start)
  G6_admr_regions$end <- G6_admr_regions$start+1
  G6_admr_regions$start <- G6_admr_regions$start*10000
  G6_admr_regions$end <- G6_admr_regions$end*10000
  G6_admr_regions$end <- as.numeric(G6_admr_regions$end)
  G6_admr_regions$start <- as.numeric(G6_admr_regions$start)
  G6_admr_regions <- select(G6_admr_regions,chr,start,end,pos)
  write.table(G6_admr_regions, "G6_10k_admr_regions.bed",sep = "\t",quote=F,row.names = F,col.names = F)
  
  
  
  
  
  
  G6_global_regions <- as.data.frame(str_split_fixed(row.names(G6_filt),"_",2))
  colnames(G6_global_regions ) <- c("chr","start")
  G6_global_regions$pos <- row.names(G6_filt)
  G6_global_regions$start <- as.numeric(G6_global_regions$start)
  G6_global_regions$end <- G6_global_regions$start+1
  G6_global_regions$start <- G6_global_regions$start*10000
  G6_global_regions$end <- G6_global_regions$end*10000
  G6_global_regions$end <- as.numeric(G6_global_regions$end)
  G6_global_regions$start <- as.numeric(G6_global_regions$start)
  G6_global_regions <- select(G6_global_regions,chr,start,end,pos)
  write.table(G6_global_regions, "G6_10k_global_regions.bed",sep = "\t",quote=F,row.names = F,col.names = F)
}





if(T){
  
  
  G6_global_H3K9me3 <- fread("g6_global_H3K9me3.bedGraph")
  G6_global_H3K9me3_count <- select(G6_global_H3K9me3,V14,V15)
  G6_global_H3K9me3_count <- aggregate(G6_global_H3K9me3_count[,2],by=list(id=G6_global_H3K9me3_count$V14),FUN=sum)
  G6_global_H3K9me3_count$hm <- 'H3K9me3'
  G6_global_H3K9me3_count$rg <- 'G6_global'
  G6_global_H3K9me3_count$V15 <- G6_global_H3K9me3_count$V15/nrow(G6_global_regions)
  
  
  G6_dmr_H3K9me3 <- fread("g6_dmr_H3K9me3.bedGraph")
  G6_dmr_H3K9me3_count <- select(G6_dmr_H3K9me3,V14,V15)
  G6_dmr_H3K9me3_count <- aggregate(G6_dmr_H3K9me3_count[,2],by=list(id=G6_dmr_H3K9me3_count$V14),FUN=sum)
  G6_dmr_H3K9me3_count$hm <- 'H3K9me3'
  G6_dmr_H3K9me3_count$rg <- 'G6_dmr'
  G6_dmr_H3K9me3_count$V15 <- G6_dmr_H3K9me3_count$V15/nrow(G6_dmr_regions)
  
  
  G6_admr_H3K9me3 <- fread("g6_admr_H3K9me3.bedGraph")
  G6_admr_H3K9me3_count <- select(G6_admr_H3K9me3,V14,V15)
  G6_admr_H3K9me3_count <- aggregate(G6_admr_H3K9me3_count[,2],by=list(id=G6_admr_H3K9me3_count$V14),FUN=sum)
  G6_admr_H3K9me3_count$hm <- 'H3K9me3'
  G6_admr_H3K9me3_count$rg <- 'G6_admr'
  G6_admr_H3K9me3_count$V15 <- G6_admr_H3K9me3_count$V15/nrow(G6_admr_regions)
  ###########################
  G6_global_H3K36me3 <- fread("g6_global_H3K36me3.bedGraph")
  G6_global_H3K36me3_count <- select(G6_global_H3K36me3,V14,V15)
  G6_global_H3K36me3_count <- aggregate(G6_global_H3K36me3_count[,2],by=list(id=G6_global_H3K36me3_count$V14),FUN=sum)
  G6_global_H3K36me3_count$hm <- 'H3K36me3'
  G6_global_H3K36me3_count$rg <- 'G6_global'
  G6_global_H3K36me3_count$V15 <- G6_global_H3K36me3_count$V15/nrow(G6_global_regions)
  
  G6_dmr_H3K36me3 <- fread("g6_dmr_H3K36me3.bedGraph")
  G6_dmr_H3K36me3_count <- select(G6_dmr_H3K36me3,V14,V15)
  G6_dmr_H3K36me3_count <- aggregate(G6_dmr_H3K36me3_count[,2],by=list(id=G6_dmr_H3K36me3_count$V14),FUN=sum)
  G6_dmr_H3K36me3_count$hm <- 'H3K36me3'
  G6_dmr_H3K36me3_count$rg <- 'G6_dmr'
  G6_dmr_H3K36me3_count$V15 <- G6_dmr_H3K36me3_count$V15/nrow(G6_dmr_regions)
  
  
  G6_admr_H3K36me3 <- fread("g6_admr_H3K36me3.bedGraph")
  G6_admr_H3K36me3_count <- select(G6_admr_H3K36me3,V14,V15)
  G6_admr_H3K36me3_count <- aggregate(G6_admr_H3K36me3_count[,2],by=list(id=G6_admr_H3K36me3_count$V14),FUN=sum)
  G6_admr_H3K36me3_count$hm <- 'H3K36me3'
  G6_admr_H3K36me3_count$rg <- 'G6_admr'
  G6_admr_H3K36me3_count$V15 <- G6_admr_H3K36me3_count$V15/nrow(G6_admr_regions)
  
  ###########################
  G6_global_H3K27me3 <- fread("g6_global_H3K27me3.bedGraph")
  G6_global_H3K27me3_count <- select(G6_global_H3K27me3,V14,V15)
  G6_global_H3K27me3_count <- aggregate(G6_global_H3K27me3_count[,2],by=list(id=G6_global_H3K27me3_count$V14),FUN=sum)
  G6_global_H3K27me3_count$hm <- 'H3K27me3'
  G6_global_H3K27me3_count$rg <- 'G6_global'
  G6_global_H3K27me3_count$V15 <- G6_global_H3K27me3_count$V15/nrow(G6_global_regions)
  
  
  G6_dmr_H3K27me3 <- fread("g6_dmr_H3K27me3.bedGraph")
  G6_dmr_H3K27me3_count <- select(G6_dmr_H3K27me3,V14,V15)
  G6_dmr_H3K27me3_count <- aggregate(G6_dmr_H3K27me3_count[,2],by=list(id=G6_dmr_H3K27me3_count$V14),FUN=sum)
  G6_dmr_H3K27me3_count$hm <- 'H3K27me3'
  G6_dmr_H3K27me3_count$rg <- 'G6_dmr'
  G6_dmr_H3K27me3_count$V15 <- G6_dmr_H3K27me3_count$V15/nrow(G6_dmr_regions)
  
  G6_admr_H3K27me3 <- fread("g6_admr_H3K27me3.bedGraph")
  G6_admr_H3K27me3_count <- select(G6_admr_H3K27me3,V14,V15)
  G6_admr_H3K27me3_count <- aggregate(G6_admr_H3K27me3_count[,2],by=list(id=G6_admr_H3K27me3_count$V14),FUN=sum)
  G6_admr_H3K27me3_count$hm <- 'H3K27me3'
  G6_admr_H3K27me3_count$rg <- 'G6_admr'
  G6_admr_H3K27me3_count$V15 <- G6_admr_H3K27me3_count$V15/nrow(G6_admr_regions)
  
  ###########################
  G6_global_H3K27ac <- fread("g6_global_H3K27ac.bedGraph")
  G6_global_H3K27ac_count <- select(G6_global_H3K27ac,V14,V15)
  G6_global_H3K27ac_count <- aggregate(G6_global_H3K27ac_count[,2],by=list(id=G6_global_H3K27ac_count$V14),FUN=sum)
  G6_global_H3K27ac_count$hm <- 'H3K27ac'
  G6_global_H3K27ac_count$rg <- 'G6_global'
  G6_global_H3K27ac_count$V15 <- G6_global_H3K27ac_count$V15/nrow(G6_global_regions)
  
  
  G6_dmr_H3K27ac <- fread("g6_dmr_H3K27ac.bedGraph")
  G6_dmr_H3K27ac_count <- select(G6_dmr_H3K27ac,V14,V15)
  G6_dmr_H3K27ac_count <- aggregate(G6_dmr_H3K27ac_count[,2],by=list(id=G6_dmr_H3K27ac_count$V14),FUN=sum)
  G6_dmr_H3K27ac_count$hm <- 'H3K27ac'
  G6_dmr_H3K27ac_count$rg <- 'G6_dmr'
  G6_dmr_H3K27ac_count$V15 <- G6_dmr_H3K27ac_count$V15/nrow(G6_dmr_regions)
  
  
  G6_admr_H3K27ac <- fread("g6_admr_H3K27ac.bedGraph")
  G6_admr_H3K27ac_count <- select(G6_admr_H3K27ac,V14,V15)
  G6_admr_H3K27ac_count <- aggregate(G6_admr_H3K27ac_count[,2],by=list(id=G6_admr_H3K27ac_count$V14),FUN=sum)
  G6_admr_H3K27ac_count$hm <- 'H3K27ac'
  G6_admr_H3K27ac_count$rg <- 'G6_admr'
  G6_admr_H3K27ac_count$V15 <- G6_admr_H3K27ac_count$V15/nrow(G6_admr_regions)
}
G6_hm_combine <- rbind(G6_global_H3K9me3_count,G6_dmr_H3K9me3_count,G6_global_H3K36me3_count,G6_dmr_H3K36me3_count,
                       G6_global_H3K27me3_count,G6_dmr_H3K27me3_count,
                       G6_global_H3K27ac_count,G6_dmr_H3K27ac_count,G6_admr_H3K27ac_count,G6_admr_H3K36me3_count,G6_admr_H3K9me3_count,G6_admr_H3K27me3_count)
colnames(G6_hm_combine) <- c("region","length","hm","region_type")
ggplot(G6_hm_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))



ggplot(G6_hm_combine)+
  geom_bar(aes(x=region_type,y=length,fill=hm),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))




G6_global_anno <- fread("g6_global_PMD.bedGraph")
G6_global_pmd <- subset(G6_global_anno,G6_global_anno$V9 == 'PMD')
G6_global_hmd <- subset(G6_global_anno,G6_global_anno$V9 == 'HMD')



G6_global_pmd_count <- select(G6_global_pmd,V4,V9,V11)
G6_global_pmd_count <- aggregate(G6_global_pmd_count[,3],by=list(id=G6_global_pmd_count$V9),FUN=sum)
G6_global_pmd_count$rg <- 'G6_global'
colnames(G6_global_pmd_count) <- c("hm","length","rg")
G6_global_pmd_count$length <- G6_global_pmd_count$length/nrow(G6_global_regions)

G6_global_hmd_count <- select(G6_global_hmd,V4,V9,V11)
G6_global_hmd_count <- aggregate(G6_global_hmd_count[,3],by=list(id=G6_global_hmd_count$V9),FUN=sum)
G6_global_hmd_count$rg <- 'G6_global'
colnames(G6_global_hmd_count) <-c("hm","length","rg")
G6_global_hmd_count$length <- G6_global_hmd_count$length/nrow(G6_global_regions)






G6_dmr_anno <- fread("g6_dmr_PMD.bedGraph")
G6_dmr_pmd <- subset(G6_dmr_anno,G6_dmr_anno$V9 == 'PMD')
G6_dmr_hmd <- subset(G6_dmr_anno,G6_dmr_anno$V9 == 'HMD')



G6_dmr_pmd_count <- select(G6_dmr_pmd,V4,V9,V11)
G6_dmr_pmd_count <- aggregate(G6_dmr_pmd_count[,3],by=list(id=G6_dmr_pmd_count$V9),FUN=sum)
G6_dmr_pmd_count$rg <- 'G6_dmr'
colnames(G6_dmr_pmd_count) <- c("hm","length","rg")
G6_dmr_pmd_count$length <- G6_dmr_pmd_count$length/nrow(G6_dmr_regions)

G6_dmr_hmd_count <- select(G6_dmr_hmd,V4,V9,V11)
G6_dmr_hmd_count <- aggregate(G6_dmr_hmd_count[,3],by=list(id=G6_dmr_hmd_count$V9),FUN=sum)
G6_dmr_hmd_count$rg <- 'G6_dmr'
colnames(G6_dmr_hmd_count) <- c("hm","length","rg")
G6_dmr_hmd_count$length <- G6_dmr_hmd_count$length/nrow(G6_dmr_regions)







G6_admr_anno <- fread("g6_admr_PMD.bedGraph")
G6_admr_pmd <- subset(G6_admr_anno,G6_admr_anno$V9 == 'PMD')
G6_admr_hmd <- subset(G6_admr_anno,G6_admr_anno$V9 == 'HMD')



G6_admr_pmd_count <- select(G6_admr_pmd,V4,V9,V11)
G6_admr_pmd_count <- aggregate(G6_admr_pmd_count[,3],by=list(id=G6_admr_pmd_count$V9),FUN=sum)
G6_admr_pmd_count$rg <- 'G6_admr'
colnames(G6_admr_pmd_count) <- c("hm","length","rg")
G6_admr_pmd_count$length <- G6_admr_pmd_count$length/nrow(G6_admr_regions)

G6_admr_hmd_count <- select(G6_admr_hmd,V4,V9,V11)
G6_admr_hmd_count <- aggregate(G6_admr_hmd_count[,3],by=list(id=G6_admr_hmd_count$V9),FUN=sum)
G6_admr_hmd_count$rg <- 'G6_admr'
colnames(G6_admr_hmd_count) <- c("hm","length","rg")
G6_admr_hmd_count$length <- G6_admr_hmd_count$length/nrow(G6_admr_regions)




G6_md_combine <- rbind(G6_admr_hmd_count,G6_admr_pmd_count,G6_dmr_hmd_count,G6_dmr_pmd_count,G6_global_hmd_count,G6_global_pmd_count)
colnames(G6_md_combine) <- c("hm","length","region_type")
ggplot(G6_md_combine)+
  geom_bar(aes(x=region_type,y=length,fill=hm),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))
ggplot(G6_md_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))





wilcox.test(scores~methods,paired=FALSE)
GA_wilcoxtest_base <- GA_test_base
for (i in 1:nrow(GA_filt)) {
  wilcox_test<-data.frame(
    ratio=unlist(GA_filt[i,]),
    group=c(rep("control",3),
            rep("oe",3))
  )
  wilcox.test_result <- wilcox.test(ratio~group,data=wilcox_test,paired=FALSE)
  GA_wilcoxtest_base[i,2] <- wilcox.test_result$p.value
  GA_wilcoxtest_base[i,3] <- mean(unlist(wilcox_test[1:3,1]))-mean(unlist(wilcox_test[4:6,1]))
}







GA_wilcoxtest_sig <- subset(GA_wilcoxtest_base,subset = GA_wilcoxtest_base$pvalue < 0.05)
GA_wilcoxtest_sig_demeth <-  subset(GA_wilcoxtest_sig,subset = GA_wilcoxtest_sig$fc > 0)
GA_wilcoxtest_sig_admeth <-  subset(GA_wilcoxtest_sig,subset = GA_wilcoxtest_sig$fc < 0)












G6baseMean = as.data.frame(colMeans(G6_filt) )
g6_pmd_regions <- G6_global_pmd$V4
G6PMDMean = as.data.frame(colMeans(G6_filt[g6_pmd_regions,]) ) 
colnames(G6PMDMean) <- "G6-PMD"
colnames(G6baseMean) <- "G6-base"





sd(x, na.rm = FALSE)
wilcox.test(scores~methods,paired=FALSE)
GA_sd_base <- GA_test_base
for (i in 1:nrow(GA_filt)) {
  GA_sd_base[i,"sd"] <- sd(GA_filt[i,], na.rm = FALSE)
}



G6_sd_base <- G6_test_base
for (i in 1:nrow(GA_filt)) {
  G6_sd_base[i,"sd"] <- sd(G6_filt[i,], na.rm = FALSE)
}





GA_sd_sort <- GA_sd_base[order(-GA_sd_base$sd),]
GA_sd_top <- GA_sd_sort[1:30000,]
GA_filt_pca <- subset(GA_filt,subset = row.names(GA_filt) %in% GA_sd_top$id)
GA_filt_pca <- as.data.frame(t(GA_filt_pca))


GA_sd_fragment <- GA_sd_sort$sd
head(GA_sd_fragment)
GA_sd_fragment <- data.frame(GA_sd_fragment)
GA_sd_res <- hist(GA_sd_fragment$GA_sd_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_sd_res$breaks),
     y = c(0, 0, GA_sd_res$counts),
     type = "l", col = "blue",
     xlab = "Standard Deviation level",
     ylab = "bin nums",
     main = "GA sd stat",
     )

abline(v=0.06227872,col="red")
abline(v=0.04396486,col="red")


GA_sd_filt_fragment <- GA_sd_top$sd
head(GA_sd_filt_fragment)
GA_sd_filt_fragment <- data.frame(GA_sd_filt_fragment)
GA_sd_filt_res <- hist(GA_sd_filt_fragment$GA_sd_filt_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, GA_sd_filt_res$breaks),
     y = c(0, 0, GA_sd_filt_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA top 30000 sd stat")

GA_group <- data.frame(Sample = rownames(GA_filt_pca), Group = rep(c("Control", "OE"), each = 3))
GA_filt_pca <- PCA(GA_filt_pca, ncp = 2, scale.unit = TRUE, graph = FALSE)
GA_filt_pca_sample <- data.frame(GA_filt_pca$ind$coord[ ,1:2])
GA_filt_pca_sample$Sample=row.names(GA_filt_pca_sample)
GA_pca_eig1 <- round(GA_filt_pca$eig[1,2], 2)
GA_pca_eig2 <- round(GA_filt_pca$eig[2,2],2 )


GA_filt_pca_sample <- merge(GA_filt_pca_sample ,GA_group,by="Sample")
head(GA_filt_pca_sample)
ggplot(data = GA_filt_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 6) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', GA_pca_eig1, '%'), y = paste('PCA2:', GA_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="GADD45A_pca")














G6_sd_sort <- G6_sd_base[order(-G6_sd_base$sd),]
G6_sd_top <- G6_sd_sort[1:30000,]
G6_filt_pca <- subset(G6_filt,subset = row.names(G6_filt) %in% G6_sd_top$id)
G6_filt_pca <- as.data.frame(t(G6_filt_pca))

G6_group <- data.frame(Sample = rownames(G6_filt_pca), Group = rep(c("Control", "OE"), each = 3))
G6_filt_pca <- PCA(G6_filt_pca, ncp = 2, scale.unit = TRUE, graph = FALSE)
G6_filt_pca_sample <- data.frame(G6_filt_pca$ind$coord[ ,1:2])
G6_filt_pca_sample$Sample=row.names(G6_filt_pca_sample)
G6_pca_eig1 <- round(G6_filt_pca$eig[1,2], 2)
G6_pca_eig2 <- round(G6_filt_pca$eig[2,2],2 )


G6_sd_fragment <- G6_sd_sort$sd
head(G6_sd_fragment)
G6_sd_fragment <- data.frame(G6_sd_fragment)
G6_sd_res <- hist(G6_sd_fragment$G6_sd_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_sd_res$breaks),
     y = c(0, 0, G6_sd_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6 sd stat")
abline(v=0.06990197,col="red")
abline(v=0.04700805,col="red")

G6_sd_filt_fragment <- G6_sd_top$sd
head(G6_sd_filt_fragment)
G6_sd_filt_fragment <- data.frame(G6_sd_filt_fragment)
G6_sd_filt_res <- hist(G6_sd_filt_fragment$G6_sd_filt_fragment, breaks = 500, plot = FALSE)
plot(x = c(0, G6_sd_filt_res$breaks),
     y = c(0, 0, G6_sd_filt_res$counts),
     type = "l", col = "blue",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6 top 30000 sd stat")


G6_filt_pca_sample <- merge(G6_filt_pca_sample ,G6_group,by="Sample")
head(G6_filt_pca_sample)
ggplot(data = G6_filt_pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 6) + 
  scale_color_manual(values = c('orange', 'purple')) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', G6_pca_eig1, '%'), y = paste('PCA2:', G6_pca_eig2, '%'), color = '') + 
  stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)+ 
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE)+
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="SNHG6_pca")




GA_r1_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[1:2], row[4:5])$p.value
})

GA_r2_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[1:2], row[5:6])$p.value
})

GA_r3_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[1:2], row[c(4,6)])$p.value
})

GA_r4_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[2:3], row[4:5])$p.value
})

GA_r5_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[2:3], row[5:6])$p.value
})

GA_r6_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[2:3], row[c(4,6)])$p.value
})

GA_r7_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[c(1,3)], row[4:5])$p.value
})

GA_r8_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[c(1,3)], row[5:6])$p.value
})

GA_r9_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[c(1,3)], row[c(4,6)])$p.value
})


GA_r10_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,6)])$p.value
})
GA_r11_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,4)])$p.value
})
GA_r12_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[c(1,6)], row[c(2,5)])$p.value
})
GA_r13_p_values <- apply(GA_filt, 1, function(row) {
  t.test(row[c(3,6)], row[c(2,5)])$p.value
})



GA_r1_test <- as.data.frame(cbind(GA_r1_p_values,row.names(GA_filt),))
GA_r2_test <- as.data.frame(cbind(GA_r2_p_values,row.names(GA_filt)))
GA_r3_test <- as.data.frame(cbind(GA_r3_p_values,row.names(GA_filt)))
GA_r4_test <- as.data.frame(cbind(GA_r4_p_values,row.names(GA_filt)))
GA_r5_test <- as.data.frame(cbind(GA_r5_p_values,row.names(GA_filt)))
GA_r6_test <- as.data.frame(cbind(GA_r6_p_values,row.names(GA_filt)))
GA_r7_test <- as.data.frame(cbind(GA_r7_p_values,row.names(GA_filt)))
GA_r8_test <- as.data.frame(cbind(GA_r8_p_values,row.names(GA_filt)))
GA_r9_test <- as.data.frame(cbind(GA_r9_p_values,row.names(GA_filt)))

GA_r10_test <- as.data.frame(cbind(GA_r10_p_values,row.names(GA_filt)))
GA_r11_test <- as.data.frame(cbind(GA_r11_p_values,row.names(GA_filt)))
GA_r12_test <- as.data.frame(cbind(GA_r12_p_values,row.names(GA_filt)))
GA_r13_test <- as.data.frame(cbind(GA_r13_p_values,row.names(GA_filt)))

GA_r1_test_sig <- subset(GA_r1_test,GA_r1_test$GA_r1_p_values < 0.05)
GA_r3_test_sig <- subset(GA_r3_test,GA_r3_test$GA_r3_p_values < 0.05)
GA_r4_test_sig <- subset(GA_r4_test,GA_r4_test$GA_r4_p_values < 0.05)
GA_r5_test_sig <- subset(GA_r5_test,GA_r5_test$GA_r5_p_values < 0.05)
GA_r6_test_sig <- subset(GA_r6_test,GA_r6_test$GA_r6_p_values < 0.05)
GA_r7_test_sig <- subset(GA_r7_test,GA_r7_test$GA_r7_p_values < 0.05)
GA_r8_test_sig <- subset(GA_r8_test,GA_r8_test$GA_r8_p_values < 0.05)
GA_r9_test_sig <- subset(GA_r9_test,GA_r9_test$GA_r9_p_values < 0.05)
GA_r10_test_sig <- subset(GA_r10_test,GA_r10_test$GA_r10_p_values < 0.05)
GA_r11_test_sig <- subset(GA_r11_test,GA_r11_test$GA_r11_p_values < 0.05)
GA_r12_test_sig <- subset(GA_r12_test,GA_r12_test$GA_r12_p_values < 0.05)
GA_r13_test_sig <- subset(GA_r13_test,GA_r13_test$GA_r13_p_values < 0.05)







G6_r1_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[1:2], row[4:5])$p.value
})
G6_r1_dt <- apply(G6_filt, 1, function(row) {
  mean(row[1:2])-mean(row[4:5])
})
G6_r2_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[1:2], row[5:6])$p.value
})
G6_r2_dt <- apply(G6_filt, 1, function(row) {
  mean(row[1:2])-mean(row[5:6])
})
G6_r3_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[1:2], row[c(4,6)])$p.value
})
G6_r3_dt <- apply(G6_filt, 1, function(row) {
  mean(row[1:2])-mean(row[c(4,6)])
})
G6_r4_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[2:3], row[4:5])$p.value
})
G6_r4_dt <- apply(G6_filt, 1, function(row) {
  mean(row[2:3])-mean(row[4:5])
})
G6_r5_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[2:3], row[5:6])$p.value
})
G6_r5_dt <- apply(G6_filt, 1, function(row) {
  mean(row[2:3])-mean(row[5:6])
})
G6_r6_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[2:3], row[c(4,6)])$p.value
})
G6_r6_dt <- apply(G6_filt, 1, function(row) {
  mean(row[2:3])-mean(row[c(4,6)])
})
G6_r7_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[c(1,3)], row[4:5])$p.value
})
G6_r7_dt <- apply(G6_filt, 1, function(row) {
  mean(row[c(1,3)])-mean(row[4:5])
})
G6_r8_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[c(1,3)], row[5:6])$p.value
})
G6_r8_dt <- apply(G6_filt, 1, function(row) {
  mean(row[c(1,3)])-mean(row[5:6])
})
G6_r9_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[c(1,3)], row[c(4,6)])$p.value
})
G6_r9_dt <- apply(G6_filt, 1, function(row) {
  mean(row[c(1,3)])-mean(row[c(4,6)])
})





G6_r10_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,6)])$p.value
})
G6_r11_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[c(1,5)], row[c(2,4)])$p.value
})
G6_r12_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[c(1,6)], row[c(2,5)])$p.value
})
G6_r13_p_values <- apply(G6_filt, 1, function(row) {
  t.test(row[c(3,6)], row[c(1,5)])$p.value
})



G6_r1_test <- as.data.frame(cbind(G6_r1_p_values,row.names(G6_filt)))
G6_r2_test <- as.data.frame(cbind(G6_r2_p_values,row.names(G6_filt)))
G6_r3_test <- as.data.frame(cbind(G6_r3_p_values,row.names(G6_filt)))
G6_r4_test <- as.data.frame(cbind(G6_r4_p_values,row.names(G6_filt)))
G6_r5_test <- as.data.frame(cbind(G6_r5_p_values,row.names(G6_filt)))
G6_r6_test <- as.data.frame(cbind(G6_r6_p_values,row.names(G6_filt)))
G6_r7_test <- as.data.frame(cbind(G6_r7_p_values,row.names(G6_filt)))
G6_r8_test <- as.data.frame(cbind(G6_r8_p_values,row.names(G6_filt)))
G6_r9_test <- as.data.frame(cbind(G6_r9_p_values,row.names(G6_filt)))
G6_r10_test <- as.data.frame(cbind(G6_r10_p_values,row.names(G6_filt)))
G6_r11_test <- as.data.frame(cbind(G6_r11_p_values,row.names(G6_filt)))
G6_r12_test <- as.data.frame(cbind(G6_r12_p_values,row.names(G6_filt)))
G6_r13_test <- as.data.frame(cbind(G6_r13_p_values,row.names(G6_filt)))



G6_r1_test_sig <- subset(G6_r1_test,G6_r1_test$G6_r1_p_values < 0.05)
G6_r2_test_sig <- subset(G6_r2_test,G6_r2_test$G6_r2_p_values < 0.05)
G6_r3_test_sig <- subset(G6_r3_test,G6_r3_test$G6_r3_p_values < 0.05)
G6_r4_test_sig <- subset(G6_r4_test,G6_r4_test$G6_r4_p_values < 0.05)
G6_r5_test_sig <- subset(G6_r5_test,G6_r5_test$G6_r5_p_values < 0.05)
G6_r6_test_sig <- subset(G6_r6_test,G6_r6_test$G6_r6_p_values < 0.05)
G6_r7_test_sig <- subset(G6_r7_test,G6_r7_test$G6_r7_p_values < 0.05)
G6_r8_test_sig <- subset(G6_r8_test,G6_r8_test$G6_r8_p_values < 0.05)
G6_r9_test_sig <- subset(G6_r9_test,G6_r9_test$G6_r9_p_values < 0.05)
G6_r10_test_sig <- subset(G6_r10_test,G6_r10_test$G6_r10_p_values < 0.05)
G6_r11_test_sig <- subset(G6_r11_test,G6_r11_test$G6_r11_p_values < 0.05)
G6_r12_test_sig <- subset(G6_r12_test,G6_r12_test$G6_r12_p_values < 0.05)
G6_r13_test_sig <- subset(G6_r13_test,G6_r13_test$G6_r13_p_values < 0.05)





GA_adjusted_p_values <- p.adjust(GA_p_values, method  = "fdr")
GA_adjusted_p_values_b <- p.adjust(GA_p_values, method  = "bonferroni")
GA_test <- cbind(GA_p_values,GA_adjusted_p_values,GA_adjusted_p_values_b,GA_dt)
GA_test <- as.data.frame(GA_test)
colnames(GA_test) <- c("pvalue","bf_adj","fdr","dt")
GA_test_significant <- subset(GA_test,subset =  pvalue< 0.05)
GA_ttest_sig_demeth <- subset(GA_test_significant , subset = dt>0)
GA_ttest_sig_admeth <- subset(GA_test_significant , subset = dt<0)


GA_test_sig <- subset(GA_test_base,subset = GA_test_base$pvalue < 0.05)
GA_test_sig_demeth <-  subset(GA_test_sig,subset = GA_test_sig$fc > 0)
GA_test_sig_admeth <-  subset(GA_test_sig,subset = GA_test_sig$fc < 0)

ga_demeth_Venn <- list(ttest = row.names(GA_ttest_sig_demeth), kwtest = GA_test_sig_demeth$id)
ga_admeth_Venn <- list(ttest = row.names(GA_ttest_sig_admeth), kwtest = GA_test_sig_admeth$id)


venn.diagram(ga_demeth_Venn, filename = 'ga_demeth_Venn.png', imagetype = 'png', 
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB'), cex = 1.5, fontfamily = 'serif')




G6baseMean[,"Sample"] <- row.names(G6baseMean)
G6baseMean[,"region"] <- "all_regions"
G6baseMean[1:3,"Sample_type"] <- "Control"
G6baseMean[4:6,"Sample_type"] <- "Overexpression"
colnames(G6baseMean) <- c("methy_level","Sample","Sample_type","region")
G6PMDMean[,"Sample"] <- row.names(G6PMDMean)
G6PMDMean[,"region"] <- "PMD_regions"
G6PMDMean[1:3,"Sample_type"] <- "Control"
G6PMDMean[4:6,"Sample_type"] <- "Overexpression"
colnames(G6PMDMean) <- c("methy_level","Sample","Sample_type","region")

G6_level <- rbind(G6baseMean,G6PMDMean)


ggplot(G6baseMean)+
  geom_boxplot(aes(x=Sample_type,y=methy_level,fill=Sample_type))+
  geom_point(aes(x=Sample_type,y=methy_level))
ggplot(G6PMDMean)+
  geom_boxplot(aes(x=Sample_type,y=methy_level,fill=Sample_type))+
  geom_point(aes(x=Sample_type,y=methy_level))

ggplot(G6baseMean)+
  geom_bar(aes(x=Sample_type,y=methy_level,fill=Sample_type),stat = "identity", position = position_dodge())+
  geom_point(aes(x=Sample_type,y=methy_level))
ggplot(G6PMDMean)+
  geom_bar(aes(x=Sample_type,y=methy_level,fill=Sample_type),stat = "identity", position = position_dodge())+
  geom_point(aes(x=Sample_type,y=methy_level))







GAbaseMean[,"region"] <- "all_regions"
GAbaseMean[1:3,"Sample_type"] <- "Control"
GAbaseMean[4:6,"Sample_type"] <- "Overexpression"
colnames(GAbaseMean) <- c("methy_level","Sample","Sample_type","region")
GAPMDMean[,"Sample"] <- row.names(GAPMDMean)
GAPMDMean[,"region"] <- "PMD_regions"
GAPMDMean[1:3,"Sample_type"] <- "Control"
GAPMDMean[4:6,"Sample_type"] <- "Overexpression"
colnames(GAPMDMean) <- c("methy_level","Sample","Sample_type","region")

GA_level <- rbind(GAbaseMean,GAPMDMean)

ggplot(GAbaseMean)+
  geom_boxplot(aes(x=Sample_type,y=methy_level,fill=Sample_type))+
  geom_point(aes(x=Sample_type,y=methy_level))
ggplot(GAPMDMean)+
  geom_boxplot(aes(x=Sample_type,y=methy_level,fill=Sample_type))+
  geom_point(aes(x=Sample_type,y=methy_level))

ggplot(GAbaseMean)+
  geom_bar(aes(x=Sample_type,y=methy_level,fill=Sample_type),stat = "identity", position = position_dodge())+
  geom_point(aes(x=Sample_type,y=methy_level))
ggplot(GAPMDMean)+
  geom_bar(aes(x=Sample_type,y=methy_level,fill=Sample_type),stat = "identity", position = position_dodge())+
  geom_point(aes(x=Sample_type,y=methy_level))



G6_diff <- as.data.frame(cbind(c(6033,5550,6169,6067,7716,5711,5162,5806),
                               c("DA","DA","DA","DA","Control","Control","Control","Control")))
colnames(G6_diff) <- c("DMR_num","group")
G6_diff$DMR_num <- as.numeric(G6_diff$DMR_num)
ggplot(G6_diff)+
  geom_bar(aes(x=group,y=DMR_num,fill=group),stat = "identity", position = position_dodge())+
  geom_point(aes(x=group,y=DMR_num))




GA_diff <- as.data.frame(cbind(c(6410,6667,6690,6907,6188,5674,6245,5913),
                               c("DA","DA","DA","DA","Control","Control","Control","Control")))
colnames(GA_diff) <- c("DMR_num","group")
GA_diff$DMR_num <- as.numeric(GA_diff$DMR_num)
ggplot(GA_diff)+
  geom_bar(aes(x=group,y=DMR_num,fill=group),stat = "identity", position = position_dodge())+
  geom_point(aes(x=group,y=DMR_num))
