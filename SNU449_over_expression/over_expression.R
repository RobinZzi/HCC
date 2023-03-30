setwd("~/projects/hcc/data/SNU449_over_expression/cpg_report")
rm(list=ls())
library(data.table)
library(R.utils)
library(tidyr)
library(pheatmap)
library(ggpointdensity)
library(dplyr)
library(viridis)

load( "methy_merge.Rdata")


save.image("methy_merge.Rdata")

fs <- list.files('./',pattern = "*.gz")

for(i in 1:length(fs)){
  data <- fread(fs[i])
  id <- str_split(paste(fs[i]),"_d")[[1]][1]
  data <- subset(data,data$V1 %in% chrs)
  colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
  data[,3] <- round(data[,2]/100)
  data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
  data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
  data[,paste(id)] <- data[,2]/(data[,2]+data[,3])
  
  if(i == 1){
    hccbase <- data[,-2]
    hccbase <- hccbase[,-2]
  }
  else{
    hccbase[,paste(id)] <- data[,4]
  }
  print(i)
}

row.names(hccbase) <- hccbase[,1]
hccbase2 <- hccbase[,-1]
hccbase3 <- na.omit(hccbase2)
pheatmap(hccbase,
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










for(i in 1:length(fs)){
  data <- fread(fs[i])
  id <- str_split(paste(fs[i]),"_d")[[1]][1]
  data <- subset(data,data$V1 %in% chrs)
  colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
  data[,3] <- round(data[,2]/10000)
  data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
  data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
  data[,paste(id)] <- data[,2]/(data[,2]+data[,3])
  
  if(i == 1){
    base_10k <- data[,-2]
    base_10k <- base_10k[,-2]
  }
  else{
    base_10k[,paste(id)] <- data[,4]
  }
  print(i)
}

base_10k_removena <- na.omit(base_10k)






ggplot(data = base, mapping = aes(x = gadd45a,
                                       y = CTNNB1)) + 
  geom_pointdensity() + #密度散点图（geom_pointdensity）
  scale_color_viridis() + 
  #geom_smooth(method = lm) +  ##省略拟合曲线
  stat_cor(method = "spearman") + 
  xlab("x") + 
  theme(axis.title.x = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  ylab("y") + 
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="serif")) +
  theme(legend.position='none')  ##去除legend

p1




base_for_dt <- base_10k
colnames(base_for_dt) <- c("chr","CTNNB1","GADD45A","GFP")
base_for_dt$ctdt <- base_for_dt$CTNNB1-base_for_dt$GFP
base_for_dt$gadt <- base_for_dt$GADD45A-base_for_dt$GFP
base_for_dt <- na.omit(base_for_dt)


ctdt <- select(base_for_dt,chr,CTNNB1,ctdt)
gadt <- select(base_for_dt,chr,GADD45A,gadt)

ctdt <- ctdt[order(ctdt$ctdt),]
gadt <- gadt[order(gadt$gadt),]




ct_subset_demeth <- ctdt[1:(0.05*nrow(ctdt)),]
ct_subset_admeth <- ctdt[(0.95*nrow(ctdt)):(nrow(ctdt)),]

ga_subset_demeth <- gadt[1:(0.05*nrow(gadt)),]
ga_subset_admeth <- gadt[(0.95*nrow(gadt)):(nrow(gadt)),]




ggplot(base_for_dt, aes(x=GADD45A, y=GFP) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_viridis() +
  stat_cor(method = "spearman") + 
  xlab("x") + 
  theme(axis.title.x = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0, 
                                    hjust = 0))+
  ylab("y") + 
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0, 
                                    hjust = 0))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="serif")) +
  theme(
    legend.position='none'
  )+
  geom_abline(slope=1, intercept=0, colour = 1, size = 0.5)



ggplot(data = base_for_dt, mapping = aes(x = CTNNB1,
                                         y = GFP)) + 
  geom_pointdensity(size=0.05) + #密度散点图（geom_pointdensity）
  scale_color_viridis() + 
  #geom_smooth(method = lm) +  ##省略拟合曲线
  stat_cor(method = "spearman") + 
  xlab("x") + 
  theme(axis.title.x = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  ylab("y") + 
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="serif")) +
  theme(legend.position='none')  +
  geom_abline(slope=1, intercept=0,se = TRUE, colour = "white", size = 0.5)


ggplot(data = base_for_dt, mapping = aes(x = GADD45A,
                                         y = GFP)) + 
  geom_pointdensity(size=0.05) + #密度散点图（geom_pointdensity）
  scale_color_viridis() + 
  #geom_smooth(method = lm) +  ##省略拟合曲线
  stat_cor(method = "spearman") + 
  xlab("x") + 
  theme(axis.title.x = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  ylab("y") + 
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="serif")) +
  theme(legend.position='none')  +
  geom_abline(slope=1, intercept=0,se = TRUE, colour = "white", size = 0.5)





ctde_graph <- as.data.frame(str_split_fixed(ct_subset_demeth$chr,"_",2))
ctde_graph$level <- as.numeric(ct_subset_demeth$CTNNB1)
colnames(ctde_graph) <- c("chr","start","level")
ctde_graph$start <- as.numeric(ctde_graph$start)
ctde_graph$end <- ctde_graph$start+1
ctde_graph$start <- ctde_graph$start*10000
ctde_graph$end <- ctde_graph$end*10000
ctde_graph$end <- as.numeric(ctde_graph$end)
ctde_graph$start <- as.numeric(ctde_graph$start)
ctde_graph <- select(ctde_graph,chr,start,end,level)
write.table(ctde_graph, "ctde.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



ctad_graph <- as.data.frame(str_split_fixed(ct_subset_admeth$chr,"_",2))
ctad_graph$level <- as.numeric(ct_subset_admeth$CTNNB1)
colnames(ctad_graph) <- c("chr","start","level")
ctad_graph$start <- as.numeric(ctad_graph$start)
ctad_graph$end <- ctad_graph$start+1
ctad_graph$start <- ctad_graph$start*10000
ctad_graph$end <- ctad_graph$end*10000
ctad_graph$end <- as.numeric(ctad_graph$end)
ctad_graph$start <- as.numeric(ctad_graph$start)
ctad_graph <- select(ctad_graph,chr,start,end,level)
write.table(ctad_graph, "ctad.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)




gade_graph <- as.data.frame(str_split_fixed(ga_subset_demeth$chr,"_",2))
gade_graph$level <- as.numeric(ga_subset_demeth$GADD45A)
colnames(gade_graph) <- c("chr","start","level")
gade_graph$start <- as.numeric(gade_graph$start)
gade_graph$end <- gade_graph$start+1
gade_graph$start <- gade_graph$start*10000
gade_graph$end <- gade_graph$end*10000
gade_graph$end <- as.numeric(gade_graph$end)
gade_graph$start <- as.numeric(gade_graph$start)
gade_graph <- select(gade_graph,chr,start,end,level)
write.table(gade_graph, "gade.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)



gaad_graph <- as.data.frame(str_split_fixed(ga_subset_admeth$chr,"_",2))
gaad_graph$level <- as.numeric(ga_subset_admeth$GADD45A)
colnames(gaad_graph) <- c("chr","start","level")
gaad_graph$start <- as.numeric(gaad_graph$start)
gaad_graph$end <- gaad_graph$start+1
gaad_graph$start <- gaad_graph$start*10000
gaad_graph$end <- gaad_graph$end*10000
gaad_graph$end <- as.numeric(gaad_graph$end)
gaad_graph$start <- as.numeric(gaad_graph$start)
gaad_graph <- select(gaad_graph,chr,start,end,level)
write.table(gaad_graph, "gaad.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)




pheatmap(base_10k_removena,
         cluster_rows = T, cluster_cols = T,
         clustering_method = "average",
         color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
         show_rownames=F,show_colnames=T,
         fontsize_col= 15,
         angle_col = 45,
         border=FALSE) 













###########ctad#################
ctad_time <- fread("ctad_time.bedGraph")
colnames(ctad_time) <- c("chr","start","end","type","chr2","start2","end2","level","length")
ctad_time <- as.data.frame(ctad_time)


ctad_time2 <- select(ctad_time,type,level,length)

ggplot()+
  geom_bar(data=ctad_time,aes(x=type,fill=type),stat = 'count')+
  theme_bw()+theme(axis.text = element_text(size = 15))

ggplot()+
  geom_bar(data=ctad_time,aes(x=type,y=length,fill=type),stat = 'summary',fun=sum)+
  theme_bw()+theme(axis.text = element_text(size = 15))



ggplot()+
  geom_bar(data=ctad_time,aes(x=type,y=level,fill=type),stat = 'summary',fun=mean)+
  theme_bw()+theme(axis.text = element_text(size = 15))
###########ctde#################
ctde_time <- fread("ctde_time.bedGraph")
colnames(ctde_time) <- c("chr","start","end","type","chr2","start2","end2","level","length")
ctde_time <- as.data.frame(ctde_time)


ggplot()+
  geom_bar(data=ctde_time,aes(x=type,fill=type),stat = 'count')+
  theme_bw()+theme(axis.text = element_text(size = 15))

ggplot()+
  geom_bar(data=ctde_time,aes(x=type,y=length,fill=type),stat = 'summary',fun=sum)+
  theme_bw()+theme(axis.text = element_text(size = 15))

ggplot()+
  geom_bar(data=ctde_time,aes(x=type,y=level,fill=type),stat = 'summary',fun=mean)+
  theme_bw()+theme(axis.text = element_text(size = 15))


###########gade#################
gade_time <- fread("gade_time.bedGraph")
colnames(gade_time) <- c("chr","start","end","type","chr2","start2","end2","level","length")
gade_time <- as.data.frame(gade_time)


ggplot()+
  geom_bar(data=gade_time,aes(x=type,fill=type),stat = 'count')+
  theme_bw()+theme(axis.text = element_text(size = 15))

ggplot()+
  geom_bar(data=gade_time,aes(x=type,y=length,fill=type),stat = 'summary',fun=sum)+
  theme_bw()+theme(axis.text = element_text(size = 15))

ggplot()+
  geom_bar(data=gade_time,aes(x=type,y=level,fill=type),stat = 'summary',fun=mean)+
  theme_bw()+theme(axis.text = element_text(size = 15))

###########gaad#################
gaad_time <- fread("gaad_time.bedGraph")
colnames(gaad_time) <- c("chr","start","end","type","chr2","start2","end2","level","length")
gaad_time <- as.data.frame(gaad_time)

ggplot()+
  geom_bar(data=gaad_time,aes(x=type,fill=type),stat = 'count')+
  theme_bw()+theme(axis.text = element_text(size = 15))

ggplot()+
  geom_bar(data=gaad_time,aes(x=type,y=length,fill=type),stat = 'summary',fun=sum)+
  theme_bw()+theme(axis.text = element_text(size = 15))

ggplot()+
  geom_bar(data=gaad_time,aes(x=type,y=level,fill=type),stat = 'summary',fun=mean)+
  theme_bw()+theme(axis.text = element_text(size = 15))






gaad_time_type <- as.data.frame(table(gaad_time$type))
gaad_time_type$level <- "hyper-meth" 
gade_time_type <- as.data.frame(table(gade_time$type))
gade_time_type$level <- "hypo-meth" 
ga_time_sum <- rbind(gaad_time_type,gade_time_type)
colnames(ga_time_sum) <- c("rep_type","count","meth_level")

gaad_time_length <- aggregate(gaad_time[,9],by=list(type=gaad_time$type),FUN=sum)
colnames(gaad_time_length) <- c("rep_type","length")
gaad_time_length$level <- "hyper-meth"
gade_time_length <- aggregate(gade_time[,9],by=list(type=gade_time$type),FUN=sum)
colnames(gade_time_length) <- c("rep_type","length")
gade_time_length$level <- "hypo-meth" 
ga_time_length_sum <- rbind(gaad_time_length,gade_time_length)
colnames(ga_time_length_sum) <- c("rep_type","length","meth_level")

gaad_time_level <- aggregate(gaad_time[,8],by=list(type=gaad_time$type),FUN=mean)
colnames(gaad_time_level) <- c("rep_type","mean_level")
gaad_time_level$level <- "hyper-meth"
gade_time_level <- aggregate(gade_time[,8],by=list(type=gade_time$type),FUN=mean)
colnames(gade_time_level) <- c("rep_type","mean_level")
gade_time_level$level <- "hypo-meth" 
ga_time_level_mean <- rbind(gaad_time_level,gade_time_level)
colnames(ga_time_level_mean) <- c("rep_type","mean_level","meth_level")

ggplot()+
  geom_bar(data=ga_time_sum,aes(x=rep_type,y=count,fill=meth_level),stat = 'identity',position = position_dodge())+
  theme_bw()+theme(axis.text = element_text(size = 15))+
  scale_fill_manual(values = c("#c7eae5",  "#af8dc3"))
ggplot()+
  geom_bar(data=ga_time_length_sum,aes(x=rep_type,y=length,fill=meth_level),stat = 'identity',position = position_dodge())+
  theme_bw()+theme(axis.text = element_text(size = 15))+
  scale_fill_manual(values = c("#c7eae5",  "#af8dc3"))
ggplot()+
  geom_bar(data=ga_time_level_mean,aes(x=rep_type,y=mean_level,fill=meth_level),stat = 'identity',position = position_dodge())+
  theme_bw()+theme(axis.text = element_text(size = 15))+
  scale_fill_manual(values = c("#c7eae5",  "#af8dc3"))



ctad_time_type <- as.data.frame(table(ctad_time$type))
ctad_time_type$level <- "hyper-meth" 
ctde_time_type <- as.data.frame(table(ctde_time$type))
ctde_time_type$level <- "hypo-meth" 
ct_time_sum <- rbind(ctad_time_type,ctde_time_type)
colnames(ct_time_sum) <- c("rep_type","count","meth_level")
ggplot()+
  geom_bar(data=ct_time_sum,aes(x=rep_type,y=count,fill=meth_level),stat = 'identity',position = position_dodge())+
  theme_bw()+theme(axis.text = element_text(size = 15))+
  scale_fill_manual(values = c("#c7eae5",  "#af8dc3"))


ggplot()+
  geom_bar(data=ga_time_sum,aes(x=rep_type,y=count,fill=meth_level),stat = 'identity',position = position_dodge())+
  theme_bw()+theme(axis.text = element_text(size = 15))+
  scale_fill_manual(values = c("#c7eae5",  "#af8dc3"))







