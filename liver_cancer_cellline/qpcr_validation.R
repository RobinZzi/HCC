
rm(list = ls())
setwd("~/projects/hcc/data/liver_cell_line/qpcr_validation")
library(pheatmap)
library(data.table)


gapdh<-as.data.frame(fread("gapdh",header = T))
row.names(gapdh) <- unlist(gapdh[,1])
gapdh <- gapdh[,-1]
gapdh_3 <- gapdh[,-4]
gapdh_3<- gapdh_3[-3,]

gapdh2<-as.data.frame(fread("gapdh2",header = T))
row.names(gapdh2) <- unlist(gapdh2[,1])
gapdh2 <- gapdh2[,-1]
gapdh2_3 <- gapdh2[,-4]
gapdh2_3<- gapdh2_3[-3,]

tubb<-as.data.frame(fread("tubb",header = T))
row.names(tubb) <- unlist(tubb[,1])
tubb <- tubb[,-1]
tubb_3 <- tubb[,-4]
tubb_3<- tubb_3[-3,]


rna_seq<-as.data.frame(fread("rna_seq",header = T))
row.names(rna_seq) <- unlist(rna_seq[,1])
rna_seq <- rna_seq[,-1]

pheatmap(gapdh_3,scale = "column",cluster_rows = F,cluster_cols = F)
pheatmap(gapdh_3,scale = "row",cluster_rows = F,cluster_cols = F)
pheatmap(gapdh_3,cluster_rows = F,cluster_cols = F)
pheatmap(gapdh_3,scale = "row",
         cluster_rows = F,cluster_cols = F,
         angle_col = 45,
         fontsize_row = 20,
         fontsize_col = 20)
pheatmap(gapdh2_3)
pheatmap(gapdh2_3,scale = "row",
         cluster_rows = F,cluster_cols = F,
         angle_col = 45,
         fontsize_row = 20,
         fontsize_col = 20)
pheatmap(tubb_3)
pheatmap(tubb_3,scale = "row",cluster_rows = F,cluster_cols = F,angle_col = 45)
pheatmap(tubb_3,scale = "row",
         cluster_rows = F,cluster_cols = F,
         angle_col = 45,
         fontsize_row = 20,
         fontsize_col = 20)
pheatmap(rna_seq)
pheatmap(rna_seq,scale = "row",cluster_rows = F,cluster_cols = F,angle_col = 45)
pheatmap(rna_seq,scale = "row",
         cluster_rows = F,cluster_cols = F,
         angle_col = 45,
         fontsize_row = 20,
         fontsize_col = 20)




g6<-counts["SNHG6",]
ga<-counts["GADD45A",]
snu_count <- select(counts,2,9,10,1)
colnames(snu_count) <- c("shGADD45A.1","shGADD45A.3","shNull.1","shNull.2")

hep_count <- select(counts,3,4,5,6,7,8)
colnames(hep_count) <- c("shGADD45A.1","shGADD45A.2","shGADD45A.3","shNull.1","shNull.2","shNull.3")




snu_result <- as.data.frame(cbind(c(62,56,24,140),
              c("shGA","shGA","shNull","shNull")))
colnames(snu_result) <- c("counts","group")
snu_result$counts <- as.numeric(snu_result$counts)
ggplot(snu_result,aes(x=group,y=counts,fill=group))+
  geom_bar(stat = 'summary',position = "dodge")+
  geom_jitter()

hep_result <- as.data.frame(cbind(c(56,2,60,174,198,160),
                                  c("shGA","shGA","shGA","shNull","shNull","shNull")))
colnames(hep_result) <- c("counts","group")
hep_result$counts <- as.numeric(hep_result$counts)
ggplot(hep_result,aes(x=group,y=counts,fill=group))+
  geom_bar(stat = 'summary',position = "dodge")+
  geom_jitter()



rpkm <- fread("combined-chrM.rpkm")
ga_rpkm <- subset(rpkm,subset = rpkm$Geneid=="GADD45A")
gapdh_rpkm <- subset(rpkm,subset = rpkm$Geneid=="GAPDH")
ga_rpkm_snu <- select(ga_rpkm,7,9,8,10,11,12)
colnames(ga_rpkm_snu) <- c("shNull.1","shNull.2","shNull.3","shGADD45A.1","shGADD45A.2","shGADD45A.3")


ga_rpkm_hep <- select(ga_rpkm,15,14,13,18,17,16)
colnames(ga_rpkm_hep) <- c("shNull.1","shNull.2","shNull.3","shGADD45A.1","shGADD45A.2","shGADD45A.3")



hep_rpkm_result <- as.data.frame(cbind(t(ga_rpkm_hep[1,]),
                                  c("shNull","shNull","shNull","shGA","shGA","shGA")))
colnames(hep_rpkm_result) <- c("rpkm","group")
hep_rpkm_result$rpkm <- as.numeric(hep_rpkm_result$rpkm)
ggplot(hep_rpkm_result,aes(x=group,y=rpkm,fill=group))+
  geom_bar(stat = 'summary',position = "dodge")+
  geom_jitter()


snu_rpkm_result <- as.data.frame(cbind(t(ga_rpkm_snu[1,]),
                                       c("shNull","shNull","shNull","shGA","shGA","shGA")))
colnames(snu_rpkm_result) <- c("rpkm","group")
snu_rpkm_result$rpkm <- as.numeric(snu_rpkm_result$rpkm)
ggplot(snu_rpkm_result,aes(x=group,y=rpkm,fill=group))+
  geom_bar(stat = 'summary',position = "dodge")+
  geom_jitter()
