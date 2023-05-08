setwd("~/projects/hcc/analysis/TCGA_related/te")

rm(list=ls())
load("te_counts.Rdata")
save.image("te_counts.Rdata")
4
4
12

anno <- fread("gdc_sample_sheet.2021-05-25.all.tsv")
anno <- anno[-1,]
liver_anno <- subset(anno ,subset = anno$`File ID` %in% bigcount$V2)

subcount <- fread("lihc.csv")

bigcount <- fread("TCGA-LIHC.csv")
bigcount <- bigcount[-1,]
colnames(bigcount) <- c("count","File ID","type")

group<- fread("methgroup.txt")
group <- group[,-1]

bigcount <- merge(bigcount,liver_anno, by = "File ID")
bigcount <- select(bigcount,count,"Sample ID")


group$sample <- gsub("\\.","-",group$sample)
colnames(group) <- c("Sample ID","group")

bigcount <- merge(bigcount,group,by="Sample ID")

total_count <- fread("TCGA-LIHC_sum.csv")
total_count <- select(total_count,)

subcount <- subcount[,-1]
subcount <- as.data.frame(subcount)
row.names(subcount) <- subcount$subfamliy
subcount <- subcount[,-1]
subcount_sum <- colSums(subcount)
subcount_sum_pm <- subcount_sum/1000000
sub_CPM

ggplot(bigcount,aes(x=group,y=count,fill=group))+
  geom_bar(stat = "identity", position = position_dodge())+theme_bw()+
  geom_point()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))






subcount <- 
