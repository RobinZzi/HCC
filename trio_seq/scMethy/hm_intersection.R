
library(ggsignif)
library(ggpubr)

save.image("methy_intersect.Rdata")

pmd_graph <- fread("big_pmd.bedGraph")
pmd_graph <- select(pmd_graph,V5,V6,V7,V9,V4)
pmd_graph <- subset(pmd_graph,subset = V9 == 'PMD')
write.table(pmd_graph, "sum_pmd.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)





 

hmd_graph <- fread("big_pmd.bedGraph")
hmd_graph <- select(hmd_graph,V5,V6,V7,V9,V4)
hmd_graph <- subset(hmd_graph,subset = V9 == 'HMD')
write.table(hmd_graph, "sum_hmd.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

















H3K27me3_pmd <- fread("H3K27me3_pmd.bedGraph")
H3K27me3_pmd_s <- select(H3K27me3_pmd,V5,V15)
H3K27me3_pmd_s$V5 <- H3K27me3_pmd$V5
H3K27me3_pmd_s <- aggregate(H3K27me3_pmd_s$V15,list(cite=H3K27me3_pmd_s$V5),sum)
colnames(H3K27me3_pmd_s) <- c("cite","length")
H3K27me3_pmd_s[,"hm_type"] <- "H3K27me3"
H3K27me3_pmd_s[,"region"] <- "pmd"
H3K27me3_hmd <- fread("H3K27me3_hmd.bedGraph")
H3K27me3_hmd_s <- select(H3K27me3_hmd,V5,V15)
H3K27me3_hmd_s$V5 <- H3K27me3_hmd$V5
H3K27me3_hmd_s <- aggregate(H3K27me3_hmd_s$V15,list(cite=H3K27me3_hmd_s$V5),sum)
colnames(H3K27me3_hmd_s) <- c("cite","length")
H3K27me3_hmd_s[,"hm_type"] <- "H3K27me3"
H3K27me3_hmd_s[,"region"] <- "hmd"


H3K27ac_pmd <- fread("H3K27ac_pmd.bedGraph")
H3K27ac_pmd_s <- select(H3K27ac_pmd,V5,V15)
H3K27ac_pmd_s$V5 <- H3K27ac_pmd$V5
H3K27ac_pmd_s <- aggregate(H3K27ac_pmd_s$V15,list(cite=H3K27ac_pmd_s$V5),sum)
colnames(H3K27ac_pmd_s) <- c("cite","length")
H3K27ac_pmd_s[,"hm_type"] <- "H3K27ac"
H3K27ac_pmd_s[,"region"] <- "pmd"
H3K27ac_hmd <- fread("H3K27ac_hmd.bedGraph")
H3K27ac_hmd_s <- select(H3K27ac_hmd,V5,V15)
H3K27ac_hmd_s$V5 <- H3K27ac_hmd$V5
H3K27ac_hmd_s <- aggregate(H3K27ac_hmd_s$V15,list(cite=H3K27ac_hmd_s$V5),sum)
colnames(H3K27ac_hmd_s) <- c("cite","length")
H3K27ac_hmd_s[,"hm_type"] <- "H3K27ac"
H3K27ac_hmd_s[,"region"] <- "hmd"


H3K9me3_pmd <- fread("H3K9me3_pmd.bedGraph")
H3K9me3_pmd_s <- select(H3K9me3_pmd,V5,V15)
H3K9me3_pmd_s$V5 <- H3K9me3_pmd$V5
H3K9me3_pmd_s <- aggregate(H3K9me3_pmd_s$V15,list(cite=H3K9me3_pmd_s$V5),sum)
colnames(H3K9me3_pmd_s) <- c("cite","length")
H3K9me3_pmd_s[,"hm_type"] <- "H3K9me3"
H3K9me3_pmd_s[,"region"] <- "pmd"
H3K9me3_hmd <- fread("H3K9me3_hmd.bedGraph")
H3K9me3_hmd_s <- select(H3K9me3_hmd,V5,V15)
H3K9me3_hmd_s$V5 <- H3K9me3_hmd$V5
H3K9me3_hmd_s <- aggregate(H3K9me3_hmd_s$V15,list(cite=H3K9me3_hmd_s$V5),sum)
colnames(H3K9me3_hmd_s) <- c("cite","length")
H3K9me3_hmd_s[,"hm_type"] <- "H3K9me3"
H3K9me3_hmd_s[,"region"] <- "hmd"



H3K4me3_pmd <- fread("H3K4me3_pmd.bedGraph")
H3K4me3_pmd_s <- select(H3K4me3_pmd,V5,V15)
H3K4me3_pmd_s$V5 <- H3K4me3_pmd$V5
H3K4me3_pmd_s <- aggregate(H3K4me3_pmd_s$V15,list(cite=H3K4me3_pmd_s$V5),sum)
colnames(H3K4me3_pmd_s) <- c("cite","length")
H3K4me3_pmd_s[,"hm_type"] <- "H3K4me3"
H3K4me3_pmd_s[,"region"] <- "pmd"
H3K4me3_hmd <- fread("H3K4me3_hmd.bedGraph")
H3K4me3_hmd_s <- select(H3K4me3_hmd,V5,V15)
H3K4me3_hmd_s$V5 <- H3K4me3_hmd$V5
H3K4me3_hmd_s <- aggregate(H3K4me3_hmd_s$V15,list(cite=H3K4me3_hmd_s$V5),sum)
colnames(H3K4me3_hmd_s) <- c("cite","length")
H3K4me3_hmd_s[,"hm_type"] <- "H3K4me3"
H3K4me3_hmd_s[,"region"] <- "hmd"


H3K36me3_pmd <- fread("H3K36me3_pmd.bedGraph")
H3K36me3_pmd_s <- select(H3K36me3_pmd,V5,V15)
H3K36me3_pmd_s$V5 <- H3K36me3_pmd$V5
H3K36me3_pmd_s <- aggregate(H3K36me3_pmd_s$V15,list(cite=H3K36me3_pmd_s$V5),sum)
colnames(H3K36me3_pmd_s) <- c("cite","length")
H3K36me3_pmd_s[,"hm_type"] <- "H3K36me3"
H3K36me3_pmd_s[,"region"] <- "pmd"
H3K36me3_hmd <- fread("H3K36me3_hmd.bedGraph")
H3K36me3_hmd_s <- select(H3K36me3_hmd,V5,V15)
H3K36me3_hmd_s$V5 <- H3K36me3_hmd$V5
H3K36me3_hmd_s <- aggregate(H3K36me3_hmd_s$V15,list(cite=H3K36me3_hmd_s$V5),sum)
colnames(H3K36me3_hmd_s) <- c("cite","length")
H3K36me3_hmd_s[,"hm_type"] <- "H3K36me3"
H3K36me3_hmd_s[,"region"] <- "hmd"


H3K9ac_pmd <- fread("H3K9ac_pmd.bedGraph")
H3K9ac_pmd_s <- select(H3K9ac_pmd,V5,V15)
H3K9ac_pmd_s$V5 <- H3K9ac_pmd$V5
H3K9ac_pmd_s <- aggregate(H3K9ac_pmd_s$V15,list(cite=H3K9ac_pmd_s$V5),sum)
colnames(H3K9ac_pmd_s) <- c("cite","length")
H3K9ac_pmd_s[,"hm_type"] <- "H3K9ac"
H3K9ac_pmd_s[,"region"] <- "pmd"
H3K9ac_hmd <- fread("H3K9ac_hmd.bedGraph")
H3K9ac_hmd_s <- select(H3K9ac_hmd,V5,V15)
H3K9ac_hmd_s$V5 <- H3K9ac_hmd$V5
H3K9ac_hmd_s <- aggregate(H3K9ac_hmd_s$V15,list(cite=H3K9ac_hmd_s$V5),sum)
colnames(H3K9ac_hmd_s) <- c("cite","length")
H3K9ac_hmd_s[,"hm_type"] <- "H3K9ac"
H3K9ac_hmd_s[,"region"] <- "hmd"



H3K27me3_sum <- rbind(H3K27me3_pmd_s,H3K27me3_hmd_s)

ggplot(H3K27me3_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()


H3K9me3_sum <- rbind(H3K9me3_pmd_s,H3K9me3_hmd_s)

ggplot(H3K9me3_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()


H3K36me3_sum <- rbind(H3K36me3_pmd_s,H3K36me3_hmd_s)

ggplot(H3K36me3_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()


H3K9ac_sum <- rbind(H3K9ac_pmd_s,H3K9ac_hmd_s)

ggplot(H3K9ac_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()


H3K27ac_sum <- rbind(H3K27ac_pmd_s,H3K27ac_hmd_s)

ggplot(H3K27ac_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()





hm_sum <- rbind(H3K27ac_sum,H3K9ac_sum,H3K36me3_sum,H3K9me3_sum,H3K27me3_sum)

ggplot(hm_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()










te_pmd <- fread("te_pmd.bedGraph")
te_pmd_s <- select(te_pmd,V5,V9,V12)
te_pmd_s_2 <- str_split_fixed(te_pmd_s$V9,":",4)
colnames(te_pmd_s_2) <- c("type","family","sub","copy")
te_pmd_s_2 <- as.data.frame(te_pmd_s_2)
te_pmd_s_3 <- select(te_pmd_s_2,type,family)
te_pmd_s <- cbind(te_pmd_s,te_pmd_s_3)
te_pmd_line <- subset(te_pmd_s,subset = type == "LINE")
te_pmd_sine <- subset(te_pmd_s,subset = type == "SINE")
te_pmd_ltr <- subset(te_pmd_s,subset = type == "LTR")
te_pmd_l1 <- subset(te_pmd_line, subset = family == "L1")
te_pmd_alu <- subset(te_pmd_sine, subset = family == "Alu")


te_pmd_l1 <- aggregate(te_pmd_l1$V12,list(cite=te_pmd_l1$V5),sum)
colnames(te_pmd_l1) <- c("cite","length")
te_pmd_l1[,"hm_type"] <- "l1"
te_pmd_l1[,"region"] <- "pmd"
te_pmd_alu <- aggregate(te_pmd_alu$V12,list(cite=te_pmd_alu$V5),sum)
colnames(te_pmd_alu) <- c("cite","length")
te_pmd_alu[,"hm_type"] <- "alu"
te_pmd_alu[,"region"] <- "pmd"

te_pmd_ltr <- aggregate(te_pmd_ltr$V12,list(cite=te_pmd_ltr$V5),sum)
colnames(te_pmd_ltr) <- c("cite","length")
te_pmd_ltr[,"hm_type"] <- "ltr"
te_pmd_ltr[,"region"] <- "pmd"




te_hmd <- fread("te_hmd.bedGraph")
te_hmd_s <- select(te_hmd,V5,V9,V12)
te_hmd_s_2 <- str_split_fixed(te_hmd_s$V9,":",4)
colnames(te_hmd_s_2) <- c("type","family","sub","copy")
te_hmd_s_2 <- as.data.frame(te_hmd_s_2)
te_hmd_s_3 <- select(te_hmd_s_2,type,family)
te_hmd_s <- cbind(te_hmd_s,te_hmd_s_3)
te_hmd_line <- subset(te_hmd_s,subset = type == "LINE")
te_hmd_sine <- subset(te_hmd_s,subset = type == "SINE")
te_hmd_ltr <- subset(te_hmd_s,subset = type == "LTR")
te_hmd_l1 <- subset(te_hmd_line, subset = family == "L1")
te_hmd_alu <- subset(te_hmd_sine, subset = family == "Alu")


te_hmd_l1 <- aggregate(te_hmd_l1$V12,list(cite=te_hmd_l1$V5),sum)
colnames(te_hmd_l1) <- c("cite","length")
te_hmd_l1[,"hm_type"] <- "l1"
te_hmd_l1[,"region"] <- "hmd"
te_hmd_alu <- aggregate(te_hmd_alu$V12,list(cite=te_hmd_alu$V5),sum)
colnames(te_hmd_alu) <- c("cite","length")
te_hmd_alu[,"hm_type"] <- "alu"
te_hmd_alu[,"region"] <- "hmd"

te_hmd_ltr <- aggregate(te_hmd_ltr$V12,list(cite=te_hmd_ltr$V5),sum)
colnames(te_hmd_ltr) <- c("cite","length")
te_hmd_ltr[,"hm_type"] <- "ltr"
te_hmd_ltr[,"region"] <- "hmd"









l1_sum <- rbind(te_hmd_l1,te_pmd_l1)
ggplot(l1_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()


alu_sum <- rbind(te_hmd_alu,te_pmd_alu)
ggplot(alu_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()


ltr_sum <- rbind(te_hmd_ltr,te_pmd_ltr)
ggplot(ltr_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot()









remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


l1_sum_removena <- l1_sum %>%
  group_by(hm_type) %>%
  mutate(length = remove_outliers(length ))
l1_sum_removena <- l1_sum_removena[complete.cases(l1_sum_removena),]
l1_sum_removena$hm_type <- gsub("l1","L1",l1_sum_removena$hm_type)
ggplot(l1_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))



alu_sum_removena <- alu_sum %>%
  group_by(hm_type) %>%
  mutate(length = remove_outliers(length ))
alu_sum_removena <- alu_sum_removena[complete.cases(alu_sum_removena),]
alu_sum_removena$hm_type <- gsub("alu","Alu",alu_sum_removena$hm_type)
ggplot(alu_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))



ltr_sum_removena <- ltr_sum %>%
  group_by(hm_type) %>%
  mutate(length = remove_outliers(length ))
ltr_sum_removena <- ltr_sum_removena[complete.cases(ltr_sum_removena),]

ggplot(ltr_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))

 
ggplot(alu_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))



ggplot(l1_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))




H3K9me3_sum <- rbind(H3K9me3_hmd_s,H3K9me3_pmd_s)
H3K9me3_sum_removena <- H3K9me3_sum %>%
  group_by(hm_type) %>%
  mutate(length = remove_outliers(length ))
H3K9me3_sum_removena <- H3K9me3_sum_removena[complete.cases(H3K9me3_sum_removena),]
ggplot(H3K9me3_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))



H3K27ac_sum <- rbind(H3K27ac_hmd_s,H3K27ac_pmd_s)
H3K27ac_sum_removena <- H3K27ac_sum %>%
  group_by(hm_type) %>%
  mutate(length = remove_outliers(length ))
H3K27ac_sum_removena <- H3K27ac_sum_removena[complete.cases(H3K27ac_sum_removena),]
ggplot(H3K27ac_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))





H3K27me3_sum <- rbind(H3K27me3_hmd_s,H3K27me3_pmd_s)
H3K27me3_sum_removena <- H3K27me3_sum %>%
  group_by(hm_type) %>%
  mutate(length = remove_outliers(length ))
H3K27me3_sum_removena <- H3K27me3_sum_removena[complete.cases(H3K27me3_sum_removena),]
ggplot(H3K27me3_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))








H3K36me3_sum <- rbind(H3K36me3_hmd_s,H3K36me3_pmd_s)
H3K36me3_sum_removena <- H3K36me3_sum %>%
  group_by(hm_type) %>%
  mutate(length = remove_outliers(length ))
H3K36me3_sum_removena <- H3K36me3_sum_removena[complete.cases(H3K36me3_sum_removena),]
ggplot(H3K36me3_sum_removena,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))





anno_sum <- rbind(H3K27ac_sum_removena,H3K36me3_sum_removena,H3K9me3_sum_removena,H3K27me3_sum_removena,ltr_sum_removena,alu_sum_removena,l1_sum_removena)
anno_sum$hm_type = factor(anno_sum$hm_type, levels=c('H3K9me3','H3K27me3','H3K27ac','H3K36me3','alu','l1','ltr'))

ggplot(anno_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))






hm_sum <- rbind(H3K27ac_sum_removena,H3K36me3_sum_removena,H3K9me3_sum_removena,H3K27me3_sum_removena)
hm_sum$hm_type = factor(hm_sum$hm_type, levels=c('H3K9me3','H3K27me3','H3K27ac','H3K36me3'))
ggplot(hm_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))+
  stat_compare_means(label = "p.signif")

he_sum <- rbind(H3K9me3_sum_removena,H3K27me3_sum_removena)
he_sum$hm_type = factor(he_sum$hm_type, levels=c('H3K9me3','H3K27me3'))
ggplot(he_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))+
  stat_compare_means(label = "p.signif")

eu_sum <- rbind(H3K36me3_sum_removena,H3K27ac_sum_removena)
eu_sum$hm_type = factor(eu_sum$hm_type, levels=c('H3K36me3','H3K27ac'))
ggplot(eu_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))+
  stat_compare_means(label = "p.signif")


te_sum <- rbind(alu_sum_removena,l1_sum_removena)
te_sum$hm_type <- gsub("alu","Alu",te_sum$hm_type)
te_sum$hm_type <- gsub("l1","L1",te_sum$hm_type)
te_sum$hm_type = factor(te_sum$hm_type, levels=c('Alu','L1'))
ggplot(te_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))+
  stat_compare_means(label = "p.signif",label.x = 1.5)

