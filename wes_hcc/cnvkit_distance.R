

rm(list = ls())

setwd("~/projects/hcc/analysis/wes_hcc")
save.image("wes_distance.Rdata")

setwd("~/projects/hcc/data/wes_hcc/hcc28/gatk_result/cnvkit_result")

hcc28_PT1 <- fread("PT1_bqsr.cnr")
hcc28_PT3 <- fread("PT3_bqsr.cnr")
hcc28_PT4 <- fread("PT4_bqsr.cnr")
hcc28_PT1_agg <- aggregate(hcc28_PT1[,6],by=list(chr=hcc28_PT1$chromosome),FUN=mean)
hcc28_PT3_agg <- aggregate(hcc28_PT3[,6],by=list(chr=hcc28_PT3$chromosome),FUN=mean)
hcc28_PT4_agg <- aggregate(hcc28_PT4[,6],by=list(chr=hcc28_PT4$chromosome),FUN=mean)
hcc28_cnv_base <- as.data.frame( cbind(hcc28_PT1_agg$log2,hcc28_PT3_agg$log2,hcc28_PT4_agg$log2))
colnames(hcc28_cnv_base) <- c("PT1","PT3","PT4")
row.names(hcc28_cnv_base) <- hcc28_PT1_agg$chr
pheatmap(hcc28_cnv_base,cluster_rows = F)

hcc28_wes_hclust <- as.data.frame(t(hcc28_cnv_base))
hcc28_wes_dist = dist(hcc28_wes_hclust, method = "euclidean")
hclust_hcc28 = hclust(hcc28_wes_dist, method = "average")
plot(hclust_hcc28)



setwd("~/projects/hcc/data/wes_hcc/hcc29/gatk_result/cnvkit_result")
hcc29_PT1 <- fread("PT1_bqsr.cnr")
hcc29_PT3 <- fread("PT3_bqsr.cnr")
hcc29_PT4 <- fread("PT4_bqsr.cnr")
hcc29_PT1_agg <- aggregate(hcc29_PT1[,6],by=list(chr=hcc29_PT1$chromosome),FUN=mean)
hcc29_PT3_agg <- aggregate(hcc29_PT3[,6],by=list(chr=hcc29_PT3$chromosome),FUN=mean)
hcc29_PT4_agg <- aggregate(hcc29_PT4[,6],by=list(chr=hcc29_PT4$chromosome),FUN=mean)
hcc29_cnv_base <- as.data.frame( cbind(hcc29_PT1_agg$log2,hcc29_PT3_agg$log2,hcc29_PT4_agg$log2))
colnames(hcc29_cnv_base) <- c("PT1","PT3","PT4")
row.names(hcc29_cnv_base) <- hcc29_PT1_agg$chr
pheatmap(hcc29_cnv_base,cluster_rows = F)

hcc29_wes_hclust <- as.data.frame(t(hcc29_cnv_base))
hcc29_wes_dist = dist(hcc29_wes_hclust, method = "euclidean")
hclust_hcc29 = hclust(hcc29_wes_dist, method = "average")
plot(hclust_hcc29)


setwd("~/projects/hcc/data/wes_hcc/hcc3/gatk_result/cnvkit_result")
hcc3_PT1 <- fread("BCPT1_bqsr.cns")
hcc3_PT2 <- fread("BCPT2_bqsr.cns")
hcc3_PT3 <- fread("BCPT3_bqsr.cns")
hcc3_PT4 <- fread("BCPT4_bqsr.cns")
hcc3_PT1_agg <- aggregate(hcc3_PT1[,5],by=list(chr=hcc3_PT1$chromosome),FUN=mean)
hcc3_PT2_agg <- aggregate(hcc3_PT2[,5],by=list(chr=hcc3_PT2$chromosome),FUN=mean)
hcc3_PT3_agg <- aggregate(hcc3_PT3[,5],by=list(chr=hcc3_PT3$chromosome),FUN=mean)
hcc3_PT4_agg <- aggregate(hcc3_PT4[,5],by=list(chr=hcc3_PT4$chromosome),FUN=mean)
hcc3_cnv_base <- as.data.frame( cbind(hcc3_PT1_agg$log2,hcc3_PT2_agg$log2,hcc3_PT3_agg$log2,hcc3_PT4_agg$log2))
colnames(hcc3_cnv_base) <- c("PT1","PT2","PT3","PT4")
row.names(hcc3_cnv_base) <- hcc3_PT1_agg$chr
pheatmap(hcc3_cnv_base,cluster_rows = F)

hcc3_wes_hclust <- as.data.frame(t(hcc3_cnv_base))
hcc3_wes_dist = dist(hcc3_wes_hclust, method = "canberra")
hclust_hcc3 = hclust(hcc3_wes_dist, method = "average")
plot(hclust_hcc3)

hcc3_PT1 <- fread("BCPT1_bqsr.cns")
hcc3_PT2 <- fread("BCPT2_bqsr.cns")
hcc3_PT3 <- fread("BCPT3_bqsr.cns")
hcc3_PT4 <- fread("BCPT4_bqsr.cns")
hcc3_PT1$wl <- hcc3_PT1$log2*hcc3_PT1$weight
hcc3_PT2$wl <- hcc3_PT2$log2*hcc3_PT2$weight
hcc3_PT3$wl <- hcc3_PT3$log2*hcc3_PT3$weight
hcc3_PT4$wl <- hcc3_PT4$log2*hcc3_PT4$weight

hcc3_PT1_agg <- aggregate(hcc3_PT1[,5],by=list(chr=hcc3_PT1$chromosome),FUN=mean)
hcc3_PT2_agg <- aggregate(hcc3_PT2[,5],by=list(chr=hcc3_PT2$chromosome),FUN=mean)
hcc3_PT3_agg <- aggregate(hcc3_PT3[,5],by=list(chr=hcc3_PT3$chromosome),FUN=mean)
hcc3_PT4_agg <- aggregate(hcc3_PT4[,5],by=list(chr=hcc3_PT4$chromosome),FUN=mean)
hcc3_cnv_base <- as.data.frame(cbind(hcc3_PT1_agg$log2,hcc3_PT2_agg$log2,hcc3_PT3_agg$log2,hcc3_PT4_agg$log2))
colnames(hcc3_cnv_base) <- c("PT1","PT2","PT3","PT4")
row.names(hcc3_cnv_base) <- hcc3_PT1_agg$chr
pheatmap(hcc3_cnv_base,cluster_rows = F,clustering_method = "centroid")

hcc3_wes_hclust <- as.data.frame(t(hcc3_cnv_base))
hcc3_wes_dist = dist(hcc3_wes_hclust, method = "euclidean")
hclust_hcc3 = hclust(hcc3_wes_dist, method = "average")
plot(hclust_hcc3)


seg <- fread("final.seg")

PT1_seg <- subset(seg, subset = V1=="BCPT1")
PT2_seg <- subset(seg, subset = V1=="BCPT2")
PT3_seg <- subset(seg, subset = V1=="BCPT3")
PT4_seg <- subset(seg, subset = V1=="BCPT4")

PT1_seg_agg <- aggregate(PT1_seg[,6],by=list(chr=PT1_seg$V2),FUN=mean)
PT2_seg_agg <- aggregate(PT2_seg[,6],by=list(chr=PT2_seg$V2),FUN=mean)
PT3_seg_agg <- aggregate(PT3_seg[,6],by=list(chr=PT3_seg$V2),FUN=mean)
PT4_seg_agg <- aggregate(PT4_seg[,6],by=list(chr=PT4_seg$V2),FUN=mean)
hcc3_cnv_base <- as.data.frame(cbind(PT1_seg_agg$V6,PT2_seg_agg$V6,PT3_seg_agg$V6,PT4_seg_agg$V6))
colnames(hcc3_cnv_base) <- c("PT1","PT2","PT3","PT4")
row.names(hcc3_cnv_base) <- PT1_seg_agg$chr
pheatmap(hcc3_cnv_base,cluster_rows = F)
