library(DSS)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/dss/HepG2")
rm(list = ls())


save.image("GA45_kd_hepg2_dss.Rdata")

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/dss/HepG2/bismark_cov")
GA45_c1 <- fread("HepG2-shNull-1.bismark.cov.gz")
GA45_c2 <- fread("HepG2-shNull-2.bismark.cov.gz")
GA45_c3 <- fread("HepG2-shNull-3.bismark.cov.gz")
GA45_e1 <- fread("HepG2-shGADD45A-1.bismark.cov.gz")
GA45_e2 <- fread("HepG2-shGADD45A-2.bismark.cov.gz")
GA45_e3 <- fread("HepG2-shGADD45A-3.bismark.cov.gz")


GA45_c1$total <- GA45_c1$V5+GA45_c1$V6
GA45_c1_con <- dplyr::select(GA45_c1,1,2,7,5)
colnames(GA45_c1_con) <- c("chr","pos","N","X")

GA45_c2$total <- GA45_c2$V5+GA45_c2$V6
GA45_c2_con <- dplyr::select(GA45_c2,1,2,7,5)
colnames(GA45_c2_con) <- c("chr","pos","N","X")

GA45_c3$total <- GA45_c3$V5+GA45_c3$V6
GA45_c3_con <- dplyr::select(GA45_c3,1,2,7,5)
colnames(GA45_c3_con) <- c("chr","pos","N","X")

GA45_e1$total <- GA45_e1$V5+GA45_e1$V6
GA45_e1_con <- dplyr::select(GA45_e1,1,2,7,5)
colnames(GA45_e1_con) <- c("chr","pos","N","X")

GA45_e2$total <- GA45_e2$V5+GA45_e2$V6
GA45_e2_con <- dplyr::select(GA45_e2,1,2,7,5)
colnames(GA45_e2_con) <- c("chr","pos","N","X")

GA45_e3$total <- GA45_e3$V5+GA45_e3$V6
GA45_e3_con <- dplyr::select(GA45_e3,1,2,7,5)
colnames(GA45_e3_con) <- c("chr","pos","N","X")

GA45_BSobj <- makeBSseqData( list(GA45_c1_con, GA45_c2_con,GA45_c3_con, GA45_e1_con, GA45_e2_con,GA45_e3_con),
                             c("C1","C2","C3", "E1", "E2","E3") )

GA45_dmlTest.sm <- DMLtest(GA45_BSobj, group1=c("C1","C2","C3"), group2=c("E1","E2","E3"), smoothing=TRUE)




GA45_dmrs <- callDMR(GA45_dmlTest.sm, p.threshold=0.01,delt=0.2,minlen = 50, minCG = 3, dis.merge = 100)

GA45_dmrs <- mutate(GA45_dmrs,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                                diff.Methy>0 ~ "Down-regulated", # 下调
                                                TRUE ~ "Unchanged"))


GA45_dmrs$region <- paste(GA45_dmrs$chr,GA45_dmrs$start,sep = "_")
table(GA45_dmrs$state)

GA45_dmls <- callDML(GA45_dmlTest.sm, p.threshold=0.01)
GA45_dmls <- mutate(GA45_dmls,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                                diff>0 ~ "Down-regulated", # 下调
                                                TRUE ~ "Unchanged"))

table(GA45_dmrs$state)
table(GA45_dmls$state)
df_color <- c("#1f76b6", "#ff7d0e")
GA45_dmrs_prop <- as.data.frame(t(table(GA45_dmrs$state)))
ggplot(GA45_dmrs_prop, aes(x = 1, y = Freq, fill = Var2))+
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

table(GA45_dmls$state)
GA45_dmls_prop <- as.data.frame(t(table(GA45_dmls$state)))
ggplot(GA45_dmls_prop, aes(x = 1, y = Freq, fill = Var2))+
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



setwd("~/projects/hcc/analysis/sh_cell_line/wgbs/dss/HepG2/dmr_anno")
write.table(GA45_dmrs, "GA45_dmrs.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
