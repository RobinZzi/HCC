setwd("~/projects/hcc/analysis/TCGA_related/methy_merge")

rm(list=ls())
save.image("late_early.RData")
load("late_early.RData")


library(ggplot2)
rank <-fread("rank.txt")
write.csv(rank,"rank.csv",sep = "/t")
if(T){
  

late_H3K9me3 <- fread("late_H3K9me3.bedGraph")
late_H3K9me3_count <- select(late_H3K9me3,V14,V15)
late_H3K9me3_count <- aggregate(late_H3K9me3_count[,2],by=list(id=late_H3K9me3_count$V14),FUN=sum)
late_H3K9me3_count$hm <- 'H3K9me3'
late_H3K9me3_count$rg <- 'late'




early_H3K9me3 <- fread("early_H3K9me3.bedGraph")
early_H3K9me3_count <- select(early_H3K9me3,V14,V15)
early_H3K9me3_count <- aggregate(early_H3K9me3_count[,2],by=list(id=early_H3K9me3_count$V14),FUN=sum)
early_H3K9me3_count$hm <- 'H3K9me3'
early_H3K9me3_count$rg <- 'early'



late_H3K36me3 <- fread("late_H3K36me3.bedGraph")
late_H3K36me3_count <- select(late_H3K36me3,V14,V15)
late_H3K36me3_count <- aggregate(late_H3K36me3_count[,2],by=list(id=late_H3K36me3_count$V14),FUN=sum)
late_H3K36me3_count$hm <- 'H3K36me3'
late_H3K36me3_count$rg <- 'late'


early_H3K36me3 <- fread("early_H3K36me3.bedGraph")
early_H3K36me3_count <- select(early_H3K36me3,V14,V15)
early_H3K36me3_count <- aggregate(early_H3K36me3_count[,2],by=list(id=early_H3K36me3_count$V14),FUN=sum)
early_H3K36me3_count$hm <- 'H3K36me3'
early_H3K36me3_count$rg <- 'early'



late_H3K4me3 <- fread("late_H3K4me3.bedGraph")
late_H3K4me3_count <- select(late_H3K4me3,V14,V15)
late_H3K4me3_count <- aggregate(late_H3K4me3_count[,2],by=list(id=late_H3K4me3_count$V14),FUN=sum)
late_H3K4me3_count$hm <- 'H3K4me3'
late_H3K4me3_count$rg <- 'late'

early_H3K4me3 <- fread("early_H3K4me3.bedGraph")
early_H3K4me3_count <- select(early_H3K4me3,V14,V15)
early_H3K4me3_count <- aggregate(early_H3K4me3_count[,2],by=list(id=early_H3K4me3_count$V14),FUN=sum)
early_H3K4me3_count$hm <- 'H3K4me3'
early_H3K4me3_count$rg <- 'early'

late_H3K27me3 <- fread("late_H3K27me3.bedGraph")
late_H3K27me3_count <- select(late_H3K27me3,V14,V15)
late_H3K27me3_count <- aggregate(late_H3K27me3_count[,2],by=list(id=late_H3K27me3_count$V14),FUN=sum)
late_H3K27me3_count$hm <- 'H3K27me3'
late_H3K27me3_count$rg <- 'late'



early_H3K27me3 <- fread("early_H3K27me3.bedGraph")
early_H3K27me3_count <- select(early_H3K27me3,V14,V15)
early_H3K27me3_count <- aggregate(early_H3K27me3_count[,2],by=list(id=early_H3K27me3_count$V14),FUN=sum)
early_H3K27me3_count$hm <- 'H3K27me3'
early_H3K27me3_count$rg <- 'early'


late_H3K27ac <- fread("late_H3K27ac.bedGraph")
late_H3K27ac_count <- select(late_H3K27ac,V14,V15)
late_H3K27ac_count <- aggregate(late_H3K27ac_count[,2],by=list(id=late_H3K27ac_count$V14),FUN=sum)
late_H3K27ac_count$hm <- 'H3K27ac'
late_H3K27ac_count$rg <- 'late'

early_H3K27ac <- fread("early_H3K27ac.bedGraph")
early_H3K27ac_count <- select(early_H3K27ac,V14,V15)
early_H3K27ac_count <- aggregate(early_H3K27ac_count[,2],by=list(id=early_H3K27ac_count$V14),FUN=sum)
early_H3K27ac_count$hm <- 'H3K27ac'
early_H3K27ac_count$rg <- 'early'

}









hm_combine <- rbind(late_H3K9me3_count,early_H3K9me3_count,late_H3K36me3_count,early_H3K36me3_count,
                 late_H3K4me3_count,early_H3K4me3_count,late_H3K27me3_count,early_H3K27me3_count,
                 late_H3K27ac_count,early_H3K27ac_count)
colnames(hm_combine) <- c("region","length","hm","region_type")
ggplot(hm_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
      axis.text=element_text(size=14,face = "bold"),
      axis.title.x=element_text(size=12),
      axis.title.y=element_text(size=13),
      legend.text=element_text(size=12),
      legend.title=element_text(size=12))



late_H3K9me3 <- fread("late_H3K9me3.bedGraph")
late_H3K9me3_count <- select(late_H3K9me3,V14,V15)
late_H3K9me3_count <- aggregate(late_H3K9me3_count[,2],by=list(id=late_H3K9me3_count$V14),FUN=sum)
late_H3K9me3_count$hm <- 'H3K9me3'
late_H3K9me3_count$rg <- 'late'




all_H3K9me3 <- fread("all_H3K9me3.bedGraph")
all_H3K9me3_count <- select(all_H3K9me3,V14,V15)
all_H3K9me3_count <- aggregate(all_H3K9me3_count[,2],by=list(id=all_H3K9me3_count$V14),FUN=sum)
all_H3K9me3_count$hm <- 'H3K9me3'
all_H3K9me3_count$rg <- 'all'



late_H3K36me3 <- fread("late_H3K36me3.bedGraph")
late_H3K36me3_count <- select(late_H3K36me3,V14,V15)
late_H3K36me3_count <- aggregate(late_H3K36me3_count[,2],by=list(id=late_H3K36me3_count$V14),FUN=sum)
late_H3K36me3_count$hm <- 'H3K36me3'
late_H3K36me3_count$rg <- 'late'


all_H3K36me3 <- fread("all_H3K36me3.bedGraph")
all_H3K36me3_count <- select(all_H3K36me3,V14,V15)
all_H3K36me3_count <- aggregate(all_H3K36me3_count[,2],by=list(id=all_H3K36me3_count$V14),FUN=sum)
all_H3K36me3_count$hm <- 'H3K36me3'
all_H3K36me3_count$rg <- 'all'



late_H3K4me3 <- fread("late_H3K4me3.bedGraph")
late_H3K4me3_count <- select(late_H3K4me3,V14,V15)
late_H3K4me3_count <- aggregate(late_H3K4me3_count[,2],by=list(id=late_H3K4me3_count$V14),FUN=sum)
late_H3K4me3_count$hm <- 'H3K4me3'
late_H3K4me3_count$rg <- 'late'

all_H3K4me3 <- fread("all_H3K4me3.bedGraph")
all_H3K4me3_count <- select(all_H3K4me3,V14,V15)
all_H3K4me3_count <- aggregate(all_H3K4me3_count[,2],by=list(id=all_H3K4me3_count$V14),FUN=sum)
all_H3K4me3_count$hm <- 'H3K4me3'
all_H3K4me3_count$rg <- 'all'

late_H3K27me3 <- fread("late_H3K27me3.bedGraph")
late_H3K27me3_count <- select(late_H3K27me3,V14,V15)
late_H3K27me3_count <- aggregate(late_H3K27me3_count[,2],by=list(id=late_H3K27me3_count$V14),FUN=sum)
late_H3K27me3_count$hm <- 'H3K27me3'
late_H3K27me3_count$rg <- 'late'



all_H3K27me3 <- fread("all_H3K27me3.bedGraph")
all_H3K27me3_count <- select(all_H3K27me3,V14,V15)
all_H3K27me3_count <- aggregate(all_H3K27me3_count[,2],by=list(id=all_H3K27me3_count$V14),FUN=sum)
all_H3K27me3_count$hm <- 'H3K27me3'
all_H3K27me3_count$rg <- 'all'


late_H3K27ac <- fread("late_H3K27ac.bedGraph")
late_H3K27ac_count <- select(late_H3K27ac,V14,V15)
late_H3K27ac_count <- aggregate(late_H3K27ac_count[,2],by=list(id=late_H3K27ac_count$V14),FUN=sum)
late_H3K27ac_count$hm <- 'H3K27ac'
late_H3K27ac_count$rg <- 'late'

all_H3K27ac <- fread("all_H3K27ac.bedGraph")
all_H3K27ac_count <- select(all_H3K27ac,V14,V15)
all_H3K27ac_count <- aggregate(all_H3K27ac_count[,2],by=list(id=all_H3K27ac_count$V14),FUN=sum)
all_H3K27ac_count$hm <- 'H3K27ac'
all_H3K27ac_count$rg <- 'all'





ha_combine <- rbind(late_H3K9me3_count,all_H3K9me3_count,late_H3K36me3_count,all_H3K36me3_count,
                    late_H3K4me3_count,all_H3K4me3_count,late_H3K27me3_count,all_H3K27me3_count,
                    late_H3K27ac_count,all_H3K27ac_count)





colnames(ha_combine) <- c("region","length","ha","region_type")
ggplot(ha_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))
  
