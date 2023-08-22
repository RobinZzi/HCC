
if(T){
  
  
  g6_global_H3K9me3 <- fread("g6_global_H3K9me3.bedGraph")
  g6_global_H3K9me3_count <- select(g6_global_H3K9me3,V14,V15)
  g6_global_H3K9me3_count <- aggregate(g6_global_H3K9me3_count[,2],by=list(id=g6_global_H3K9me3_count$V14),FUN=sum)
  g6_global_H3K9me3_count$hm <- 'H3K9me3'
  g6_global_H3K9me3_count$rg <- 'g6_global'
  g6_global_H3K9me3_count$V15 <- g6_global_H3K9me3_count$V15/nrow(g6_global_regions)
  
  
  g6_dmr_H3K9me3 <- fread("g6_dmr_H3K9me3.bedGraph")
  g6_dmr_H3K9me3_count <- select(g6_dmr_H3K9me3,V14,V15)
  g6_dmr_H3K9me3_count <- aggregate(g6_dmr_H3K9me3_count[,2],by=list(id=g6_dmr_H3K9me3_count$V14),FUN=sum)
  g6_dmr_H3K9me3_count$hm <- 'H3K9me3'
  g6_dmr_H3K9me3_count$rg <- 'g6_dmr'
  g6_dmr_H3K9me3_count$V15 <- g6_dmr_H3K9me3_count$V15/nrow(g6_dmr_regions)
  
  
  g6_admr_H3K9me3 <- fread("g6_admr_H3K9me3.bedGraph")
  g6_admr_H3K9me3_count <- select(g6_admr_H3K9me3,V14,V15)
  g6_admr_H3K9me3_count <- aggregate(g6_admr_H3K9me3_count[,2],by=list(id=g6_admr_H3K9me3_count$V14),FUN=sum)
  g6_admr_H3K9me3_count$hm <- 'H3K9me3'
  g6_admr_H3K9me3_count$rg <- 'g6_admr'
  g6_admr_H3K9me3_count$V15 <- g6_admr_H3K9me3_count$V15/nrow(g6_admr_regions)
  ###########################
  g6_global_H3K36me3 <- fread("g6_global_H3K36me3.bedGraph")
  g6_global_H3K36me3_count <- select(g6_global_H3K36me3,V14,V15)
  g6_global_H3K36me3_count <- aggregate(g6_global_H3K36me3_count[,2],by=list(id=g6_global_H3K36me3_count$V14),FUN=sum)
  g6_global_H3K36me3_count$hm <- 'H3K36me3'
  g6_global_H3K36me3_count$rg <- 'g6_global'
  g6_global_H3K36me3_count$V15 <- g6_global_H3K36me3_count$V15/nrow(g6_global_regions)
  
  g6_dmr_H3K36me3 <- fread("g6_dmr_H3K36me3.bedGraph")
  g6_dmr_H3K36me3_count <- select(g6_dmr_H3K36me3,V14,V15)
  g6_dmr_H3K36me3_count <- aggregate(g6_dmr_H3K36me3_count[,2],by=list(id=g6_dmr_H3K36me3_count$V14),FUN=sum)
  g6_dmr_H3K36me3_count$hm <- 'H3K36me3'
  g6_dmr_H3K36me3_count$rg <- 'g6_dmr'
  g6_dmr_H3K36me3_count$V15 <- g6_dmr_H3K36me3_count$V15/nrow(g6_dmr_regions)
  

  g6_admr_H3K36me3 <- fread("g6_admr_H3K36me3.bedGraph")
  g6_admr_H3K36me3_count <- select(g6_admr_H3K36me3,V14,V15)
  g6_admr_H3K36me3_count <- aggregate(g6_admr_H3K36me3_count[,2],by=list(id=g6_admr_H3K36me3_count$V14),FUN=sum)
  g6_admr_H3K36me3_count$hm <- 'H3K36me3'
  g6_admr_H3K36me3_count$rg <- 'g6_admr'
  g6_admr_H3K36me3_count$V15 <- g6_admr_H3K36me3_count$V15/nrow(g6_admr_regions)
    
  ###########################
  g6_global_H3K27me3 <- fread("g6_global_H3K27me3.bedGraph")
  g6_global_H3K27me3_count <- select(g6_global_H3K27me3,V14,V15)
  g6_global_H3K27me3_count <- aggregate(g6_global_H3K27me3_count[,2],by=list(id=g6_global_H3K27me3_count$V14),FUN=sum)
  g6_global_H3K27me3_count$hm <- 'H3K27me3'
  g6_global_H3K27me3_count$rg <- 'g6_global'
  g6_global_H3K27me3_count$V15 <- g6_global_H3K27me3_count$V15/nrow(g6_global_regions)
  
  
  g6_dmr_H3K27me3 <- fread("g6_dmr_H3K27me3.bedGraph")
  g6_dmr_H3K27me3_count <- select(g6_dmr_H3K27me3,V14,V15)
  g6_dmr_H3K27me3_count <- aggregate(g6_dmr_H3K27me3_count[,2],by=list(id=g6_dmr_H3K27me3_count$V14),FUN=sum)
  g6_dmr_H3K27me3_count$hm <- 'H3K27me3'
  g6_dmr_H3K27me3_count$rg <- 'g6_dmr'
  g6_dmr_H3K27me3_count$V15 <- g6_dmr_H3K27me3_count$V15/nrow(g6_dmr_regions)
  
  g6_admr_H3K27me3 <- fread("g6_admr_H3K27me3.bedGraph")
  g6_admr_H3K27me3_count <- select(g6_admr_H3K27me3,V14,V15)
  g6_admr_H3K27me3_count <- aggregate(g6_admr_H3K27me3_count[,2],by=list(id=g6_admr_H3K27me3_count$V14),FUN=sum)
  g6_admr_H3K27me3_count$hm <- 'H3K27me3'
  g6_admr_H3K27me3_count$rg <- 'g6_admr'
  g6_admr_H3K27me3_count$V15 <- g6_admr_H3K27me3_count$V15/nrow(g6_admr_regions)
  
  ###########################
  g6_global_H3K27ac <- fread("g6_global_H3K27ac.bedGraph")
  g6_global_H3K27ac_count <- select(g6_global_H3K27ac,V14,V15)
  g6_global_H3K27ac_count <- aggregate(g6_global_H3K27ac_count[,2],by=list(id=g6_global_H3K27ac_count$V14),FUN=sum)
  g6_global_H3K27ac_count$hm <- 'H3K27ac'
  g6_global_H3K27ac_count$rg <- 'g6_global'
  g6_global_H3K27ac_count$V15 <- g6_global_H3K27ac_count$V15/nrow(g6_global_regions)
  
  
  g6_dmr_H3K27ac <- fread("g6_dmr_H3K27ac.bedGraph")
  g6_dmr_H3K27ac_count <- select(g6_dmr_H3K27ac,V14,V15)
  g6_dmr_H3K27ac_count <- aggregate(g6_dmr_H3K27ac_count[,2],by=list(id=g6_dmr_H3K27ac_count$V14),FUN=sum)
  g6_dmr_H3K27ac_count$hm <- 'H3K27ac'
  g6_dmr_H3K27ac_count$rg <- 'g6_dmr'
  g6_dmr_H3K27ac_count$V15 <- g6_dmr_H3K27ac_count$V15/nrow(g6_dmr_regions)
  
  
  g6_admr_H3K27ac <- fread("g6_admr_H3K27ac.bedGraph")
  g6_admr_H3K27ac_count <- select(g6_admr_H3K27ac,V14,V15)
  g6_admr_H3K27ac_count <- aggregate(g6_admr_H3K27ac_count[,2],by=list(id=g6_admr_H3K27ac_count$V14),FUN=sum)
  g6_admr_H3K27ac_count$hm <- 'H3K27ac'
  g6_admr_H3K27ac_count$rg <- 'g6_admr'
  g6_admr_H3K27ac_count$V15 <- g6_admr_H3K27ac_count$V15/nrow(g6_admr_regions)
}
g6_hm_combine <- rbind(g6_global_H3K9me3_count,g6_dmr_H3K9me3_count,g6_global_H3K36me3_count,g6_dmr_H3K36me3_count,
                       g6_global_H3K27me3_count,g6_dmr_H3K27me3_count,
                       g6_global_H3K27ac_count,g6_dmr_H3K27ac_count,g6_admr_H3K27ac_count,g6_admr_H3K36me3_count,g6_admr_H3K9me3_count,g6_admr_H3K27me3_count)
colnames(g6_hm_combine) <- c("region","length","hm","region_type")
ggplot(g6_hm_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))


g6_hm_combine <- rbind(g6_admr_H3K9me3_count,g6_dmr_H3K9me3_count,g6_admr_H3K36me3_count,g6_dmr_H3K36me3_count,
                       g6_admr_H3K27me3_count,g6_dmr_H3K27me3_count,
                       g6_admr_H3K27ac_count,g6_dmr_H3K27ac_count)
colnames(g6_hm_combine) <- c("region","length","hm","region_type")
ggplot(g6_hm_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))





if(T){
  
  
  ga_admr_H3K9me3 <- fread("ga_admr_H3K9me3.bedGraph")
  ga_admr_H3K9me3_count <- select(ga_admr_H3K9me3,V14,V15)
  ga_admr_H3K9me3_count <- aggregate(ga_admr_H3K9me3_count[,2],by=list(id=ga_admr_H3K9me3_count$V14),FUN=sum)
  ga_admr_H3K9me3_count$hm <- 'H3K9me3'
  ga_admr_H3K9me3_count$rg <- 'ga_admr'
  
  
  
  
  ga_dmr_H3K9me3 <- fread("ga_dmr_H3K9me3.bedGraph")
  ga_dmr_H3K9me3_count <- select(ga_dmr_H3K9me3,V14,V15)
  ga_dmr_H3K9me3_count <- aggregate(ga_dmr_H3K9me3_count[,2],by=list(id=ga_dmr_H3K9me3_count$V14),FUN=sum)
  ga_dmr_H3K9me3_count$hm <- 'H3K9me3'
  ga_dmr_H3K9me3_count$rg <- 'ga_dmr'
  
  
  
  ga_admr_H3K36me3 <- fread("ga_admr_H3K36me3.bedGraph")
  ga_admr_H3K36me3_count <- select(ga_admr_H3K36me3,V14,V15)
  ga_admr_H3K36me3_count <- aggregate(ga_admr_H3K36me3_count[,2],by=list(id=ga_admr_H3K36me3_count$V14),FUN=sum)
  ga_admr_H3K36me3_count$hm <- 'H3K36me3'
  ga_admr_H3K36me3_count$rg <- 'ga_admr'
  
  
  ga_dmr_H3K36me3 <- fread("ga_dmr_H3K36me3.bedGraph")
  ga_dmr_H3K36me3_count <- select(ga_dmr_H3K36me3,V14,V15)
  ga_dmr_H3K36me3_count <- aggregate(ga_dmr_H3K36me3_count[,2],by=list(id=ga_dmr_H3K36me3_count$V14),FUN=sum)
  ga_dmr_H3K36me3_count$hm <- 'H3K36me3'
  ga_dmr_H3K36me3_count$rg <- 'ga_dmr'
  
  
  
  
  ga_admr_H3K27me3 <- fread("ga_admr_H3K27me3.bedGraph")
  ga_admr_H3K27me3_count <- select(ga_admr_H3K27me3,V14,V15)
  ga_admr_H3K27me3_count <- aggregate(ga_admr_H3K27me3_count[,2],by=list(id=ga_admr_H3K27me3_count$V14),FUN=sum)
  ga_admr_H3K27me3_count$hm <- 'H3K27me3'
  ga_admr_H3K27me3_count$rg <- 'ga_admr'
  
  
  
  ga_dmr_H3K27me3 <- fread("ga_dmr_H3K27me3.bedGraph")
  ga_dmr_H3K27me3_count <- select(ga_dmr_H3K27me3,V14,V15)
  ga_dmr_H3K27me3_count <- aggregate(ga_dmr_H3K27me3_count[,2],by=list(id=ga_dmr_H3K27me3_count$V14),FUN=sum)
  ga_dmr_H3K27me3_count$hm <- 'H3K27me3'
  ga_dmr_H3K27me3_count$rg <- 'ga_dmr'
  
  
  ga_admr_H3K27ac <- fread("ga_admr_H3K27ac.bedGraph")
  ga_admr_H3K27ac_count <- select(ga_admr_H3K27ac,V14,V15)
  ga_admr_H3K27ac_count <- aggregate(ga_admr_H3K27ac_count[,2],by=list(id=ga_admr_H3K27ac_count$V14),FUN=sum)
  ga_admr_H3K27ac_count$hm <- 'H3K27ac'
  ga_admr_H3K27ac_count$rg <- 'ga_admr'
  
  ga_dmr_H3K27ac <- fread("ga_dmr_H3K27ac.bedGraph")
  ga_dmr_H3K27ac_count <- select(ga_dmr_H3K27ac,V14,V15)
  ga_dmr_H3K27ac_count <- aggregate(ga_dmr_H3K27ac_count[,2],by=list(id=ga_dmr_H3K27ac_count$V14),FUN=sum)
  ga_dmr_H3K27ac_count$hm <- 'H3K27ac'
  ga_dmr_H3K27ac_count$rg <- 'ga_dmr'
  
}
ga_hm_combine <- rbind(ga_admr_H3K9me3_count,ga_dmr_H3K9me3_count,ga_admr_H3K36me3_count,ga_dmr_H3K36me3_count,
                       ga_admr_H3K27me3_count,ga_dmr_H3K27me3_count,
                       ga_admr_H3K27ac_count,ga_dmr_H3K27ac_count)
colnames(ga_hm_combine) <- c("region","length","hm","region_type")
ggplot(ga_hm_combine)+
  geom_bar(aes(x=hm,y=length,fill=region_type),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))











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


if(T){
  g6_global_anno <- fread("g6_global_PMD.bedGraph")
  g6_global_pmd <- subset(g6_global_anno,g6_global_anno$V9 == 'PMD')
  g6_global_hmd <- subset(g6_global_anno,g6_global_anno$V9 == 'HMD')
  
  
  
  g6_global_pmd_count <- select(g6_global_pmd,V4,V9,V11)
  g6_global_pmd_count <- aggregate(g6_global_pmd_count[,3],by=list(id=g6_global_pmd_count$V9),FUN=sum)
  g6_global_pmd_count$rg <- 'g6_global'
  colnames(g6_global_pmd_count) <- c("hm","length","rg")
  g6_global_pmd_count$length <- g6_global_pmd_count$length/nrow(g6_global_regions)
  
  g6_global_hmd_count <- select(g6_global_hmd,V4,V9,V11)
  g6_global_hmd_count <- aggregate(g6_global_hmd_count[,3],by=list(id=g6_global_hmd_count$V9),FUN=sum)
  g6_global_hmd_count$rg <- 'g6_global'
  colnames(g6_global_hmd_count) <-c("hm","length","rg")
  g6_global_hmd_count$length <- g6_global_hmd_count$length/nrow(g6_global_regions)
  
  
  
  
  
  
  g6_dmr_anno <- fread("g6_dmr_PMD.bedGraph")
  g6_dmr_pmd <- subset(g6_dmr_anno,g6_dmr_anno$V9 == 'PMD')
  g6_dmr_hmd <- subset(g6_dmr_anno,g6_dmr_anno$V9 == 'HMD')
  
  
  
  g6_dmr_pmd_count <- select(g6_dmr_pmd,V4,V9,V11)
  g6_dmr_pmd_count <- aggregate(g6_dmr_pmd_count[,3],by=list(id=g6_dmr_pmd_count$V9),FUN=sum)
  g6_dmr_pmd_count$rg <- 'g6_dmr'
  colnames(g6_dmr_pmd_count) <- c("hm","length","rg")
  g6_dmr_pmd_count$length <- g6_dmr_pmd_count$length/nrow(g6_dmr_regions)
  
  g6_dmr_hmd_count <- select(g6_dmr_hmd,V4,V9,V11)
  g6_dmr_hmd_count <- aggregate(g6_dmr_hmd_count[,3],by=list(id=g6_dmr_hmd_count$V9),FUN=sum)
  g6_dmr_hmd_count$rg <- 'g6_dmr'
  colnames(g6_dmr_hmd_count) <- c("hm","length","rg")
  g6_dmr_hmd_count$length <- g6_dmr_hmd_count$length/nrow(g6_dmr_regions)
  
  
  
  
  
  
  
  g6_admr_anno <- fread("g6_admr_PMD.bedGraph")
  g6_admr_pmd <- subset(g6_admr_anno,g6_admr_anno$V9 == 'PMD')
  g6_admr_hmd <- subset(g6_admr_anno,g6_admr_anno$V9 == 'HMD')
  
  
  
  g6_admr_pmd_count <- select(g6_admr_pmd,V4,V9,V11)
  g6_admr_pmd_count <- aggregate(g6_admr_pmd_count[,3],by=list(id=g6_admr_pmd_count$V9),FUN=sum)
  g6_admr_pmd_count$rg <- 'g6_admr'
  colnames(g6_admr_pmd_count) <- c("hm","length","rg")
  g6_admr_pmd_count$length <- g6_admr_pmd_count$length/nrow(g6_admr_regions)
  
  g6_admr_hmd_count <- select(g6_admr_hmd,V4,V9,V11)
  g6_admr_hmd_count <- aggregate(g6_admr_hmd_count[,3],by=list(id=g6_admr_hmd_count$V9),FUN=sum)
  g6_admr_hmd_count$rg <- 'g6_admr'
  colnames(g6_admr_hmd_count) <- c("hm","length","rg")
  g6_admr_hmd_count$length <- g6_admr_hmd_count$length/nrow(g6_admr_regions)
}



g6_md_combine <- rbind(g6_admr_hmd_count,g6_admr_pmd_count,g6_dmr_hmd_count,g6_dmr_pmd_count,g6_global_hmd_count,g6_global_pmd_count)
colnames(g6_md_combine) <- c("hm","length","region_type")
ggplot(g6_md_combine)+
  geom_bar(aes(x=region_type,y=length,fill=hm),stat = "identity", position = position_dodge())+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))


