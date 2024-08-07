library(DSS)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs-20240611/dss/HepG2")
rm(list = ls())




setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240611/meth_ex")
SNHG6_c1 <- fread("SNHG6_V1_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
SNHG6_c2 <- fread("SNHG6_V2_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
SNHG6_e1 <- fread("SNHG6_SH2-1_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240614/meth_ex")
SNHG6_e2 <- fread("SNHG6_E2_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")

save.image("dss.Rdata")

SNHG6_c1$total <- SNHG6_c1$V5+SNHG6_c1$V6
SNHG6_c1_con <- dplyr::select(SNHG6_c1,1,2,7,5)
colnames(SNHG6_c1_con) <- c("chr","pos","N","X")

SNHG6_c2$total <- SNHG6_c2$V5+SNHG6_c2$V6
SNHG6_c2_con <- dplyr::select(SNHG6_c2,1,2,7,5)
colnames(SNHG6_c2_con) <- c("chr","pos","N","X")

SNHG6_e1$total <- SNHG6_e1$V5+SNHG6_e1$V6
SNHG6_e1_con <- dplyr::select(SNHG6_e1,1,2,7,5)
colnames(SNHG6_e1_con) <- c("chr","pos","N","X")

SNHG6_e2$total <- SNHG6_e2$V5+SNHG6_e2$V6
SNHG6_e2_con <- dplyr::select(SNHG6_e2,1,2,7,5)
colnames(SNHG6_e2_con) <- c("chr","pos","N","X")

SNG6_BSobj <- makeBSseqData( list(SNHG6_c1_con, SNHG6_c2_con, SNHG6_e1_con, SNHG6_e2_con),
                        c("C1","C2", "E1", "E2") )

SNG6_dmlTest.sm <- DMLtest(SNG6_BSobj, group1=c("C1", "C2"), group2=c("E1", "E2"), smoothing=TRUE)


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

setwd("~/projects/hcc/analysis/sh_cell_line/wgbs_20240611/methy_mtx/dmr_anno")
write.table(SNG6_dmrs, "SNG6_dmrs.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)

