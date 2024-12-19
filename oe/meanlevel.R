data1 <-  rbind(GAHMDMean$`GA-HMD`,GAPMDMean$`GA-PMD`,GAbaseMean$`GA-global`)





G6meansum <- cbind(G6HMDMean,G6PMDMean,G6baseMean)
GAmeansum$sample <- c(rep("V",3),rep("OE",3))
G6meansum$sample <- c(rep("V",3),rep("OE",3))

colnames(GAmeansum) <- c("HMD_Methy","PMD_Methy","Global_Methy","Group") 
colnames(G6meansum) <- c("HMD_Methy","PMD_Methy","Global_Methy","Group") 


data1 <- as.data.frame(cbind(c(G6meansum$Global_Methy,G6meansum$HMD_Methy,G6meansum$PMD_Methy),
                             c(rep("Global_Mean",6),rep("PMD_Mean",6),rep("HMD_Mean",6)),
                             c(rep(c(rep("V",3),rep("OE",3)),3))))
colnames(data1) <- c("Methy_level","Region","Group") 

data1$Methy_level <- as.numeric(data1$Methy_level)

ggplot(data1,aes(Group,Methy_level,fill=Group))+
  theme_bw()+
  geom_bar(stat = "identity",position =position_dodge(0.1))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=0,angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=00,color = 'black'),
        axis.title.y = element_text(size=00,color = 'black'),
        axis.title.x = element_text(size=0,color = 'black'))+
  stat_summary(fun = mean, geom = "bar",position =position_dodge(0.1))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",width = 0.2,position = position_dodge(0.1))+
  geom_jitter()+
  facet_wrap(~ Region)+ylim(0,1)
















data4 <- as.data.frame(cbind(c(GAmeansum$Global_Methy,GAmeansum$HMD_Methy,GAmeansum$PMD_Methy),
                             c(rep("Global_Mean",6),rep("PMD_Mean",6),rep("HMD_Mean",6)),
                             c(rep(c(rep("V",3),rep("OE",3)),3))))
colnames(data4) <- c("Methy_level","Region","Group") 

data4$Methy_level <- as.numeric(data4$Methy_level)

ggplot(data2,aes(Group,Methy_level,fill=Group))+
  theme_bw()+
  geom_bar(stat = "identity",position =position_dodge(0.1))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=0,angle = 45,hjust = 1,color = 'black'),
        axis.text.y = element_text(size=00,color = 'black'),
        axis.title.y = element_text(size=00,color = 'black'),
        axis.title.x = element_text(size=0,color = 'black'))+
  stat_summary(fun = mean, geom = "bar",position =position_dodge(0.1))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",width = 0.2,position = position_dodge(0.1))+
  geom_jitter()+
  facet_wrap(~ Region)+ylim(0,1)
