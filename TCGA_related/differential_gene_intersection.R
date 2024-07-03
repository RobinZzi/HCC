setwd("~/projects/hcc/analysis/TCGA_related/rna_intersect")
rm(list=ls())

save.image("diff_intersect.RData")
load("~/projects/hcc/analysis/TCGA_related/diff_intersect.RData")


TCGA_diff <- as.data.frame(fread("DEG_edgeR_hypoall_vs_hyper.txt"))
row.names(TCGA_diff) <- TCGA_diff[,1]
TCGA_diff <- TCGA_diff[,-1]

hcc3_demeth_down <- as.data.frame(fread("hcc3_demeth_down.txt"))
row.names(hcc3_demeth_down) <- hcc3_demeth_down[,1]
hcc3_demeth_down <- hcc3_demeth_down[,-1]
hcc3_demeth_down_gene <- hcc3_demeth_down$gene

hcc3_demeth_up <- as.data.frame(fread("hcc3_demeth_up.txt"))
row.names(hcc3_demeth_up) <- hcc3_demeth_up[,1]
hcc3_demeth_up <- hcc3_demeth_up[,-1]
hcc3_demeth_up_gene <- hcc3_demeth_up$gene

hcc28_demeth_up <- as.data.frame(fread("hcc28_demeth_up.txt"))
row.names(hcc28_demeth_up) <- hcc28_demeth_up[,1]
hcc28_demeth_up <- hcc28_demeth_up[,-1]
hcc28_demeth_up_gene <- hcc28_demeth_up$gene

hcc28_demeth_down <- as.data.frame(fread("hcc28_demeth_down.txt"))
row.names(hcc28_demeth_down) <- hcc28_demeth_down[,1]
hcc28_demeth_down <- hcc28_demeth_down[,-1]
hcc28_demeth_down_gene <- hcc28_demeth_down$gene


hcc29_demeth_up <- as.data.frame(fread("hcc29_demeth_up.txt"))
row.names(hcc29_demeth_up) <- hcc29_demeth_up[,1]
hcc29_demeth_up <- hcc29_demeth_up[,-1]
hcc29_demeth_up_gene <- hcc29_demeth_up$gene

hcc29_demeth_down <- as.data.frame(fread("hcc29_demeth_down.txt"))
row.names(hcc29_demeth_down) <- hcc29_demeth_down[,1]
hcc29_demeth_down <- hcc29_demeth_down[,-1]
hcc29_demeth_down_gene <- hcc29_demeth_down$gene


DEG_edgeR_hypoup <- subset(TCGA_diff, subset = TCGA_diff$change == 'up')
DEG_edgeR_hypodown <- subset(TCGA_diff, subset = TCGA_diff$change == 'down')


TCGA_hypoup_gene <- DEG_edgeR_hypoup$gene
TCGA_hypodown_gene <- DEG_edgeR_hypodown$gene



sc_common_up <- intersect(hcc3_demeth_up_gene,hcc28_demeth_up_gene)
sc_common_up <- intersect(sc_common_up,hcc29_demeth_up_gene)

sc_common_down <- intersect(hcc3_demeth_down_gene,hcc28_demeth_down_gene)
sc_common_down <- intersect(sc_common_down,hcc29_demeth_down_gene)


sc_tcga_common_up <- intersect(sc_common_up,TCGA_hypoup_gene)
sc_tcga_common_down <- intersect(sc_common_down,TCGA_hypodown_gene)


sc28_tcga_common_down <- intersect(hcc28_demeth_down_gene,TCGA_hypodown_gene)
sc29_tcga_common_down <- intersect(hcc29_demeth_down_gene,TCGA_hypodown_gene)
sc3_tcga_common_down <- intersect(hcc3_demeth_down_gene,TCGA_hypodown_gene)

sc28_tcga_common_up <- intersect(hcc28_demeth_up_gene,TCGA_hypoup_gene)
sc29_tcga_common_up <- intersect(hcc29_demeth_up_gene,TCGA_hypoup_gene)
sc3_tcga_common_up <- intersect(hcc3_demeth_up_gene,TCGA_hypoup_gene)










up_venn <- list(hcc3_demeth_up_gene,hcc28_demeth_up_gene,hcc29_demeth_up_gene,TCGA_hypoup_gene)
up_venn2 <- list(hcc3_demeth_up_gene,hcc28_demeth_up_gene,hcc29_demeth_up_gene)

venn.diagram(up_venn, filename = 'up.png', imagetype = 'png', ,category.names = c("hcc3","hcc28","hcc29","TCGA"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#61C8FF'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#61C8FF'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#61C8FF'), cex = 1.5, fontfamily = 'serif')
venn.diagram(up_venn2, filename = 'up2.png', imagetype = 'png', ,category.names = c("hcc3","hcc28","hcc29"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')



down_venn <- list(hcc3_demeth_down_gene,hcc28_demeth_down_gene,hcc29_demeth_down_gene,TCGA_hypodown_gene)

venn.diagram(down_venn, filename = 'down.png', imagetype = 'png', ,category.names = c("hcc3","hcc28","hcc29","TCGA"),
             fill = c('#4D157D', '#84C7DB', '#C8D948','#61C8FF'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#61C8FF'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#61C8FF'), cex = 1.5, fontfamily = 'serif')

