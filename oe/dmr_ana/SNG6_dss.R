SNG6_c1 <- fread("SNHG6-C1.bismark.cov.gz")
SNG6_c2 <- fread("SNHG6-C2.bismark.cov.gz")
SNG6_c3 <- fread("SNHG6-C3.bismark.cov.gz")
SNG6_e1 <- fread("SNHG6-OE1.bismark.cov.gz")
SNG6_e2 <- fread("SNHG6-OE2.bismark.cov.gz")
SNG6_e3 <- fread("SNHG6-OE3.bismark.cov.gz")

SNG6_c1$total <- SNG6_c1$V5+SNG6_c1$V6
SNG6_c1_con <- dplyr::select(SNG6_c1,1,2,7,5)
colnames(SNG6_c1_con) <- c("chr","pos","N","X")

SNG6_c2$total <- SNG6_c2$V5+SNG6_c2$V6
SNG6_c2_con <- dplyr::select(SNG6_c2,1,2,7,5)
colnames(SNG6_c2_con) <- c("chr","pos","N","X")

SNG6_c3$total <- SNG6_c3$V5+SNG6_c3$V6
SNG6_c3_con <- dplyr::select(SNG6_c3,1,2,7,5)
colnames(SNG6_c3_con) <- c("chr","pos","N","X")

SNG6_e1$total <- SNG6_e1$V5+SNG6_e1$V6
SNG6_e1_con <- dplyr::select(SNG6_e1,1,2,7,5)
colnames(SNG6_e1_con) <- c("chr","pos","N","X")

SNG6_e2$total <- SNG6_e2$V5+SNG6_e2$V6
SNG6_e2_con <- dplyr::select(SNG6_e2,1,2,7,5)
colnames(SNG6_e2_con) <- c("chr","pos","N","X")

SNG6_e3$total <- SNG6_e3$V5+SNG6_e3$V6
SNG6_e3_con <- dplyr::select(SNG6_e3,1,2,7,5)
colnames(SNG6_e3_con) <- c("chr","pos","N","X")
SNG6_BSobj <- makeBSseqData( list(SNG6_c1_con, SNG6_c2_con,SNG6_c3_con, SNG6_e1_con, SNG6_e2_con,SNG6_e3_con),
                             c("C1","C2","C3", "E1", "E2","E3") )

SNG6_dmlTest.sm <- DMLtest(SNG6_BSobj, group1=c("C1","C2","C3"), group2=c("E1","E2","E3"), smoothing=TRUE)



SNG6_dmrs <- callDMR(SNG6_dmlTest.sm, p.threshold=0.01,delt=0.2,minlen = 50, minCG = 3, dis.merge = 100)

SNG6_dmrs <- mutate(SNG6_dmrs,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                           diff.Methy>0 ~ "Down-regulated", # 下调
                                           TRUE ~ "Unchanged"))
table(SNG6_dmrs$state)

SNG6_dmls <- callDML(SNG6_dmlTest.sm, p.threshold=0.01)
SNG6_dmls <- mutate(SNG6_dmls,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                                diff>0 ~ "Down-regulated", # 下调
                                                TRUE ~ "Unchanged"))


table(SNG6_dmls$state)
showOneDMR(SNG6_dmrs['1001',], SNG6_BSobj)
df_color <- c("#1f76b6", "#ff7d0e")
SNG6_dmrs_prop <- as.data.frame(t(table(SNG6_dmrs$state)))
ggplot(SNG6_dmrs_prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))

table(SNG6_dmls$state)
SNG6_dmls_prop <- as.data.frame(t(table(SNG6_dmls$state)))
ggplot(SNG6_dmls_prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))

setwd("~/projects/hcc/analysis/oe/dss/SNG6_oe")
write.table(SNG6_dmrs, "SNG6_dmrs.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

