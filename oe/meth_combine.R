setwd("~/projects/hcc/analysis/oe/merge/merged_level")
rm(list=ls())

GA_C1 <- fread("GA_C1_SUM.txt")
GA_C2 <- fread("GA_C2_SUM.txt")
GA_C3 <- fread("GA_C3_SUM.txt")
GA_E1 <- fread("GA_OE1_SUM.txt")
GA_E2 <- fread("GA_OE2_SUM.txt")
GA_E3 <- fread("GA_OE3_SUM.txt")


GA <- cbind(GA_C1,GA_C2[,5],GA_C3[,5],GA_E1[,5],GA_E2[,5],GA_E3[,5])
GA_mtx <- as.data.frame(GA[,5:10])
row.names(GA_mtx) <- GA$id
GA_mtx <- na.omit(GA_mtx)
GA_geatmap <- pheatmap(GA_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,treeheight_row = 0,
         angle_col = 45,annotation_row = pmd_row_anno,
         border=FALSE)
GAbaseMean = as.data.frame(colMeans(GA_mtx) ) 


G6_C1 <- fread("G6_C1_SUM.txt")
G6_C2 <- fread("G6_C2_SUM.txt")
G6_C3 <- fread("G6_C3_SUM.txt")
G6_E1 <- fread("G6_OE1_SUM.txt")
G6_E2 <- fread("G6_OE2_SUM.txt")
G6_E3 <- fread("G6_OE3_SUM.txt")

G6 <- cbind(G6_C1,G6_C2[,5],G6_C3[,5],G6_E1[,5],G6_E2[,5],G6_E3[,5])
G6_mtx <- as.data.frame(G6[,5:10])
row.names(G6_mtx) <- G6$id
G6_mtx <- na.omit(G6_mtx)
pheatmap(G6_mtx,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "canberra",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,treeheight_row = 0,
         angle_col = 45,annotation_row = pmd_row_anno,
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
G6PMDMean = as.data.frame(colMeans(G6_mtx[pmd_regions,]) ) 
colnames(G6PMDMean) <- "G6-PMD"
colnames(G6baseMean) <- "G6-base"
colnames(GAbaseMean) <- "GA-base"



GA_C1$sumcount <- GA_C1$meth+GA_C1$demeth
GA_C2$sumcount <- GA_C1$meth+GA_C2$demeth
GA_C3$sumcount <- GA_C1$meth+GA_C3$demeth
GA_E1$sumcount <- GA_C1$meth+GA_E1$demeth
GA_E2$sumcount <- GA_C1$meth+GA_E2$demeth
GA_E3$sumcount <- GA_C1$meth+GA_E3$demeth
G6_C1$sumcount <- G6_C1$meth+G6_C1$demeth
G6_C2$sumcount <- G6_C1$meth+G6_C2$demeth
G6_C3$sumcount <- G6_C1$meth+G6_C3$demeth
G6_E1$sumcount <- G6_C1$meth+G6_E1$demeth
G6_E2$sumcount <- G6_C1$meth+G6_E2$demeth
G6_E3$sumcount <- G6_C1$meth+G6_E3$demeth

GA_C1$sample <- "GA_C1"
GA_C2$sample <- "GA_C2"
GA_C3$sample <- "GA_C3"
GA_E1$sample <- "GA_E1"
GA_E2$sample <- "GA_E2"
GA_E3$sample <- "GA_E3"
G6_C1$sample <- "G6_C1"
G6_C2$sample <- "G6_C2"
G6_C3$sample <- "G6_C3"
G6_E1$sample <- "G6_E1"
G6_E2$sample <- "G6_E2"
G6_E3$sample <- "G6_E3"


GA_C1_count <- select(GA_C1,id,sumcount,sample)
GA_C2_count <- select(GA_C2,id,sumcount,sample)
GA_C3_count <- select(GA_C3,id,sumcount,sample)
GA_E1_count <- select(GA_E1,id,sumcount,sample)
GA_E2_count <- select(GA_E2,id,sumcount,sample)
GA_E3_count <- select(GA_E3,id,sumcount,sample)
G6_C1_count <- select(G6_C1,id,sumcount,sample)
G6_C2_count <- select(G6_C2,id,sumcount,sample)
G6_C3_count <- select(G6_C3,id,sumcount,sample)
G6_E1_count <- select(G6_E1,id,sumcount,sample)
G6_E2_count <- select(G6_E2,id,sumcount,sample)
G6_E3_count <- select(G6_E3,id,sumcount,sample)



GA_count <- rbind(GA_C1_count,GA_C2_count,GA_C3_count,GA_E1_count,GA_E2_count,GA_E3_count)
G6_count <- rbind(G6_C1_count,G6_C2_count,G6_C3_count,G6_E1_count,G6_E2_count,G6_E3_count)

ggplot(GA_count,aes(x=id,y=sumcount,fill=sample))+
  geom_bar(stat = "identity", position = position_dodge(0.5),width = 0.4)








length_cutoff <- 5000
breaks_num <- 5000

GA_fragment <- GA_count$sumcount[GA_count$sumcount <= length_cutoff]
head(GA_fragment)
GA_fragment <- data.frame(GA_fragment)
GA_res <- hist(GA_fragment$GA_fragment, breaks = breaks_num, plot = FALSE)
plot(x = c(0, GA_res$breaks),
     y = c(0, 0, GA_res$counts),
     type = "l", col = "red",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "GA Sample reads stat")

ggplot(GA_fragment, aes(x = GA_fragment)) + 
  geom_histogram(aes(y = ..count..),
                 colour = 1, fill = "white") +theme_bw()

G6_fragment <- G6_count$sumcount[G6_count$sumcount <= length_cutoff]
head(G6_fragment)
G6_fragment <- data.frame(G6_fragment)
G6_res <- hist(G6_fragment$G6_fragment, breaks = breaks_num, plot = FALSE)
plot(x = c(0, G6_res$breaks),
     y = c(0, 0, G6_res$counts),
     type = "l", col = "red",
     xlab = "Sumed_count",
     ylab = "bin nums",
     main = "G6 Sample reads stat")

ggplot(G6_fragment, aes(x = G6_fragment)) + 
  geom_histogram(aes(y = ..count..),
                 colour = 1, fill = "white") +theme_bw()











ga_order_row = GA_geatmap$tree_row$order
ga_order_col = GA_geatmap$tree_col$order
sort_ga_all_data = data.frame(GA_mtx[ga_order_row, ga_order_col])
ga_demeth <- sort_ga_all_data[13100:15590,]
pheatmap(ga_demeth,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,treeheight_row = 0,
         angle_col = 45,annotation_row = pmd_row_anno,
         border=FALSE)
ga_dmr_regions <- as.data.frame(str_split_fixed(row.names(ga_demeth),"_",2))
colnames(ga_dmr_regions) <- c("chr","start")
ga_dmr_regions$pos <- row.names(ga_demeth)
ga_dmr_regions$start <- as.numeric(ga_dmr_regions$start)
ga_dmr_regions$end <- ga_dmr_regions$start+1
ga_dmr_regions$start <- ga_dmr_regions$start*100000
ga_dmr_regions$end <- ga_dmr_regions$end*100000
ga_dmr_regions$end <- as.numeric(ga_dmr_regions$end)
ga_dmr_regions$start <- as.numeric(ga_dmr_regions$start)
ga_dmr_regions <- select(ga_dmr_regions,chr,start,end,pos)
write.table(ga_dmr_regions, "ga_dmr_regions.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)




dmr_H3K27ac <- fread("dmr_H3K27ac.bedGraph")
dmr_H3K27ac_count <- select(dmr_H3K27ac,V14,V15)
dmr_H3K27ac_count <- aggregate(dmr_H3K27ac_count[,2],by=list(id=dmr_H3K27ac_count$V14),FUN=sum)
dmr_H3K27ac_count$hm <- 'H3K27ac'
dmr_H3K27ac_count$rg <- 'dmr'

dmr_H3K36me3 <- fread("dmr_H3K36me3.bedGraph")
dmr_H3K36me3_count <- select(dmr_H3K36me3,V14,V15)
dmr_H3K36me3_count <- aggregate(dmr_H3K36me3_count[,2],by=list(id=dmr_H3K36me3_count$V14),FUN=sum)
dmr_H3K36me3_count$hm <- 'H3K36me3'
dmr_H3K36me3_count$rg <- 'dmr'

dmr_H3K9me3 <- fread("dmr_H3K9me3.bedGraph")
dmr_H3K9me3_count <- select(dmr_H3K9me3,V14,V15)
dmr_H3K9me3_count <- aggregate(dmr_H3K9me3_count[,2],by=list(id=dmr_H3K9me3_count$V14),FUN=sum)
dmr_H3K9me3_count$hm <- 'H3K9me3'
dmr_H3K9me3_count$rg <- 'dmr'


dmr_H3K27me3 <- fread("dmr_H3K27me3.bedGraph")
dmr_H3K27me3_count <- select(dmr_H3K27me3,V14,V15)
dmr_H3K27me3_count <- aggregate(dmr_H3K27me3_count[,2],by=list(id=dmr_H3K27me3_count$V14),FUN=sum)
dmr_H3K27me3_count$hm <- 'H3K27me3'
dmr_H3K27me3_count$rg <- 'dmr'

hm_combine <- rbind(dmr_H3K27me3_count,dmr_H3K9me3_count,dmr_H3K36me3_count,dmr_H3K27ac_count)
colnames(hm_combine) <- c("region","length","hm","region_type")
ggplot(hm_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))
