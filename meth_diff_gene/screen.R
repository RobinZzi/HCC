setwd("~/projects/hcc/analysis/meth_diff_gene/screen")
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
pot <- function(x) {
  2^x
}

save.image("screen.Rdata")
meth_delta <- fread("methy_change.txt")

pt1_diff <- fread("hcc2_pt1_diff.txt")

pt2_diff <- fread("hcc2_pt2_diff.txt")

common_up <- fread("pt1_pt2_up.txt")

common_down <- fread("pt1_pt2_down.txt")

pt1_commonup <- subset(pt1_diff, subset=symbol %in% common_up$x)

pt1_commondown <- subset(pt1_diff, subset=symbol %in% common_down$x)

pt2_commonup <- subset(pt2_diff, subset=symbol %in% common_up$x)

pt2_commondown <- subset(pt2_diff, subset=symbol %in% common_down$x)

pt1_commonup <- dplyr::select(pt1_commonup,1,3)
pt2_commonup <- dplyr::select(pt2_commonup,1,3)
colnames(pt1_commonup) <- c("gene","pt1")
colnames(pt2_commonup) <- c("gene","pt2")
common_up_fc <- merge(pt1_commonup,pt2_commonup,by='gene')
common_up_fc$rna_change <- 'up'
  
pt1_commondown <- dplyr::select(pt1_commondown,1,3)
pt2_commondown <- dplyr::select(pt2_commondown,1,3)
colnames(pt1_commondown) <- c("gene","pt1")
colnames(pt2_commondown) <- c("gene","pt2")
common_down_fc <- merge(pt1_commondown,pt2_commondown,by='gene') 
common_down_fc$rna_change <- 'down'
sum <- rbind(common_up_fc,common_down_fc)

meth_delta_filt <- subset(meth_delta,subset=gene %in% sum_up$gene)
sum <- subset(sum,subset=gene %in% meth_delta_filt$gene)
meth_delta_filt <- dplyr::select(meth_delta_filt,2,6)
sum <- merge(sum,meth_delta_filt,by='gene')

sum <- dplyr::mutate(sum,avg_fc = log2((pot(pt1)+pot(pt2))/2))
sum <- dplyr::mutate(sum,label = case_when(gene %in% c("SNHG6","GADD45A") ~ gene))
sum <- dplyr::mutate(sum,color = case_when(gene %in% c("SNHG6","GADD45A") ~ 'red'))
sum <- dplyr::select(sum,1,5,6,7,8)
colnames(sum) <- c("gene","avg_meth_dt","log_avg_rna_FC","label","color")






ggplot(sum,aes(x=avg_meth_dt,y=log_avg_rna_FC))+
  geom_point(size=0.5,color="grey")+
  geom_point(data = sum[sum$label != "", ], color = "red")+
  theme_bw()+
  geom_text_repel(aes(label=label, color = color), size =5,hjust=3,vjust=-3)+
  theme(legend.position = "none")
  
