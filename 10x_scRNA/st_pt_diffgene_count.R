data <- as.data.frame(cbind(c(410,32,10,4,22,508,0,381,0,0),
                            rep(c("HCC1","HCC7","HCC6","HCC5","HCC3"),2),
                            c(rep(c("Up"),5),rep(c("Down"),5))))
colnames(data) <- c("Num","Patient","Type")                      

data$Num <- as.numeric(data$Num)

data$Type <- factor(data$Type,levels = c("Up","Down"))

library(ggplot2)

library(ggpubr)
ggplot(data,aes(Patient,Num,fill=Type))+
  geom_bar(width=0.5,stat = 'identity',position = position_dodge(0.5))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank())
