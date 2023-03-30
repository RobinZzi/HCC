setwd("~/projects/hcc/data/SNU449_over_expression/cpg_report") 

rm(list=ls())
library(data.table)
library(R.utils)
library(tidyr)
library(pheatmap)
library(edgeR)


load("methy_merge.Rdata")


save.image("methy_merge.Rdata")

    data <- fread("SNU449-GADD45A_remove_dup.CpG_report.txt.gz")
    data <- subset(data,data$V1 %in% chrs)
    colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
    data[,3] <- round(data[,2]/100000)
    data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
    data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
    data[,'gadd45a'] <- data[,2]/(data[,2]+data[,3])
    base <- data[,-2]
    base <- base[,-2]
    
    data <- fread("SNU449-CTNNB1_remove_dup.CpG_report.txt.gz")
    data <- subset(data,data$V1 %in% chrs)
    colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
    data[,3] <- round(data[,2]/100000)
    data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
    data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
    data[,'CTNNB1'] <- data[,2]/(data[,2]+data[,3])
    base[,'CTNNB1'] <- data[,4]
    
    data <- fread("SNU449-GFP_remove_dup.CpG_report.txt.gz")
    data <- subset(data,data$V1 %in% chrs)
    colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
    data[,3] <- round(data[,2]/100000)
    data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
    data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
    data[,'GFP'] <- data[,2]/(data[,2]+data[,3])
    base[,'GFP'] <- data[,4]

    base_heatmap <- base
    row.names(base_heatmap) <- base[,1]
    base_heatmap <- base_heatmap[,-1]
    base_heatmap <- as.numeric(base_heatmap)
    base_heatmap2 <- as.data.frame(lapply(base_heatmap,as.numeric))    
    row.names(base_heatmap2) <- row.names(base_heatmap)
    
    base_heatmap3 <- na.omit(base_heatmap2) 
    colnames(base_heatmap3) <- c("GADD45A","CTNNB1","GFP")
    
    pheatmap(base_heatmap3,
             cluster_rows = T, cluster_cols = T,
             clustering_method = "average",
             color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
             show_rownames=F,show_colnames=T,
             fontsize_col= 15,
             angle_col = 45,
             border=FALSE) 
    baseMean = colMeans(base_heatmap3)    
    
    
    
    
    
    ga_graph <- as.data.frame(str_split_fixed(row.names(base_heatmap3),"_",2))
    ga_graph$count <- as.numeric(base_heatmap3$GADD45A)
    colnames(ga_graph) <- c("chr","start","count")
    ga_graph$start <- as.numeric(ga_graph$start)
    ga_graph$end <- ga_graph$start+1
    ga_graph$start <- ga_graph$start*100000
    ga_graph$end <- ga_graph$end*100000
    ga_graph$end <- as.numeric(ga_graph$end)
    ga_graph$start <- as.numeric(ga_graph$start)
    ga_graph <- select(ga_graph,chr,start,end,count)
    write.table(ga_graph, "ga.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
    
    
    
    
    gfp_graph <- as.data.frame(str_split_fixed(row.names(base_heatmap3),"_",2))
    gfp_graph$count <- as.numeric(base_heatmap3$GFP)
    colnames(gfp_graph) <- c("chr","start","count")
    gfp_graph$start <- as.numeric(gfp_graph$start)
    gfp_graph$end <- gfp_graph$start+1
    gfp_graph$start <- gfp_graph$start*100000
    gfp_graph$end <- gfp_graph$end*100000
    gfp_graph$end <- as.numeric(gfp_graph$end)
    gfp_graph$start <- as.numeric(gfp_graph$start)
    gfp_graph <- select(gfp_graph,chr,start,end,count)
    write.table(gfp_graph, "gfp.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
    