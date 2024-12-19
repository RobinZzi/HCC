library(ggplot2)
save.image("rtpcr.Rdata")

oe_g6_datas <- c("0.000188128","0.001243826","0.019077289")
oe_g6_datas <- as.numeric(oe_g6_datas)
oe_g6_samples <- c("HepG2_wt","HepG2_SNHG6-DOX(-)","HepG2_SNHG6-DOX(+)") 
oe_g6 <- as.data.frame(cbind(oe_g6_samples,oe_g6_datas))

colnames(oe_g6) <- c("samples","Relative Expression Level")

oe_g6$`Relative Expression Level` <- as.numeric(oe_g6$`Relative Expression Level`)
oe_g6$samples <- factor(oe_g6$samples,levels = c("HepG2_wt","HepG2_SNHG6-DOX(-)","HepG2_SNHG6-DOX(+)"))
ggplot(data=oe_g6)+
  geom_bar(aes(x=samples,y=`Relative Expression Level`,fill=samples),stat='summary',fun=mean)+
  theme_bw()+
  scale_fill_manual(values =c('#b11a2b','#4a74a4','#4aa44d'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=18,face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=16,face="bold",color = 'black'),
        axis.title.y = element_text(size=16,face="bold",color = 'black'),
        axis.title.x = element_text(size=0,face="bold",color = 'black'))


oe_ga_datas <- c("0.006407601","0.013221228","0.990342872")
oe_ga_datas <- as.numeric(oe_ga_datas)
oe_ga_samples <- c("HepG2_wt","HepG2_GADD45A-DOX(-)","HepG2_GADD45A-DOX(+)") 
oe_ga <- as.data.frame(cbind(oe_ga_samples,oe_ga_datas))

colnames(oe_ga) <- c("samples","Relative Expression Level")

oe_ga$`Relative Expression Level` <- as.numeric(oe_ga$`Relative Expression Level`)
oe_ga$samples <- factor(oe_ga$samples,levels = c("HepG2_wt","HepG2_GADD45A-DOX(-)","HepG2_GADD45A-DOX(+)"))
ggplot(data=oe_ga)+
  geom_bar(aes(x=samples,y=`Relative Expression Level`,fill=samples),stat='summary',fun=mean)+
  theme_bw()+NoLegend()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=20,face="bold",color = 'black'),
        axis.title.y = element_text(size=30,face="bold",color = 'black'),
        axis.title.x = element_text(size=0,face="bold",color = 'black'))

ggplot(data=oe_ga)+
  geom_bar(aes(x=samples,y=`Relative Expression Level`,fill=samples),stat='summary',fun=mean)+
  theme_bw()+
  scale_fill_manual(values =c('#b11a2b','#4a74a4','#4aa44d'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=18,face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=16,face="bold",color = 'black'),
        axis.title.y = element_text(size=16,face="bold",color = 'black'),
        axis.title.x = element_text(size=0,face="bold",color = 'black'))


library(ggplot2)

sh_g6_datas <- c("0.009427579","0.006116311","0.003488884")
sh_g6_datas <- as.numeric(sh_g6_datas)
sh_g6_samples <- c("Huh7_control","Huh7_SNHG6_KD1","Huh7_SNHG6_KD2") 
sh_g6 <- as.data.frame(cbind(sh_g6_samples,sh_g6_datas))

colnames(sh_g6) <- c("samples","Relative Expression Level")

sh_g6$`Relative Expression Level` <- as.numeric(sh_g6$`Relative Expression Level`)
sh_g6$samples <- factor(sh_g6$samples,levels = c("Huh7_control","Huh7_SNHG6_KD1","Huh7_SNHG6_KD2") )
ggplot(data=sh_g6)+
  geom_bar(aes(x=samples,y=`Relative Expression Level`,fill=samples),stat='summary',fun=mean)+
  theme_bw()+NoLegend()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=16,face="bold",color = 'black'),
        axis.title.y = element_text(size=30,face="bold",color = 'black'),
        axis.title.x = element_text(size=0,face="bold",color = 'black'))


sh_ga_datas <- c("0.009427579","0.006116311","0.003488884")
sh_ga_datas <- as.numeric(sh_ga_datas)
sh_ga_samples <- c("293T_control","293T_GADD45A_KD1","293T_GADD45A_KD2") 
sh_ga <- as.data.frame(cbind(sh_ga_samples,sh_ga_datas))

colnames(sh_ga) <- c("samples","Relative Expression Level")

sh_ga$`Relative Expression Level` <- as.numeric(sh_ga$`Relative Expression Level`)
sh_ga$samples <- factor(sh_ga$samples,levels = c("293T_control","293T_GADD45A_KD1","293T_GADD45A_KD2"))
ggplot(data=sh_ga)+
  geom_bar(aes(x=samples,y=`Relative Expression Level`,fill=samples),stat='summary',fun=mean)+
  theme_bw()+NoLegend()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=20,face="bold",color = 'black'),
        axis.title.y = element_text(size=30,face="bold",color = 'black'),
        axis.title.x = element_text(size=0,face="bold",color = 'black'))



library(data.table)
library(ggplot2)
library(ggpubr)



KD_data <- fread("qpcr-sh")
GA_KD <- KD_data[1:6,]
G6_KD <- KD_data[7:12,]

GA_KD$group <- c("V","V","CRi1","CRi1","CRi2","CRi2")
colnames(GA_KD) <- c("sample","level","group")
GA_KD$group <- factor(GA_KD$group,levels=c("V","CRi1","CRi2"))
ggplot(data=GA_KD,aes(x=group,y=level,color=group))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  geom_jitter(width = 0.1,size=2)+
  scale_color_manual(values =c('#b11a2b','#4a74a4','#4aa44d'))+
  stat_compare_means(comparisons = list(c("V","CRi1")),
                     method = "t.test",label = "p.signif",
                     label.y =1 )+
  stat_compare_means(comparisons = list(c("V","CRi2")),
                     method = "t.test",label = "p.signif",
                     label.y =1.1 )+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=20,face="bold",color = 'black'),
        axis.title.y = element_text(size=30,face="bold",color = 'black'),
        axis.title.x = element_text(size=0,face="bold",color = 'black'))


G6_KD$group <- c("V","V","SH1","SH1","SH2","SH2")
colnames(G6_KD) <- c("sample","level","group")
G6_KD$group <- factor(G6_KD$group,levels=c("V","SH1","SH2"))
ggplot(data=G6_KD,aes(x=group,y=level,fill=group))+
  theme_bw()+
  geom_bar(stat = "identity",position =position_dodge(0.1))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4','#4aa44d'))+
  stat_compare_means(comparisons = list(c("V","SH1")),
                                                method = "t.test", 
                                                label.y =0.01 )+
  stat_compare_means(comparisons = list(c("V","SH2")),
                     method = "t.test", 
                     label.y =0.011 )+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=20,color = 'black'),
        axis.title.y = element_text(size=30,color = 'black'),
        axis.title.x = element_text(size=0,color = 'black'))+
  stat_summary(fun = mean, geom = "bar",position =position_dodge(0.1))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",width = 0.2,position = position_dodge(0.1))+
  geom_jitter()



ggplot(data=GA_KD,aes(x=group,y=level,fill=group))+
  theme_bw()+
  geom_bar(stat = "identity",position =position_dodge(0.1))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4','#4aa44d'))+
  stat_compare_means(comparisons = list(c("V","CRi1")),
                     method = "t.test",
                     label.y = 1 )+
  stat_compare_means(comparisons = list(c("V","CRi2")),
                     method = "t.test",
                     label.y = 1.1 )+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=20,color = 'black'),
        axis.title.y = element_text(size=30,color = 'black'),
        axis.title.x = element_text(size=0,color = 'black'))+
  stat_summary(fun = mean, geom = "bar",position =position_dodge(0.1))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",width = 0.2,position = position_dodge(0.1))+
  geom_jitter()



G6_OE <- oe_g6
colnames(G6_OE) <- colnames(G6_KD[,1:2])

ggplot(data=G6_OE,aes(x=sample,y=level,fill=sample))+
  theme_bw()+
  geom_bar(stat = "identity",position =position_dodge(0.1))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4','#4aa44d'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=20,color = 'black'),
        axis.title.y = element_text(size=30,color = 'black'),
        axis.title.x = element_text(size=0,color = 'black'))


GA_OE <- oe_ga
colnames(GA_OE) <- colnames(GA_KD[,1:2])

ggplot(data=GA_OE,aes(x=sample,y=level,fill=sample))+
  theme_bw()+
  geom_bar(stat = "identity",position =position_dodge(0.1))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4','#4aa44d'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=22,angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=20,color = 'black'),
        axis.title.y = element_text(size=30,color = 'black'),
        axis.title.x = element_text(size=0,color = 'black'))
