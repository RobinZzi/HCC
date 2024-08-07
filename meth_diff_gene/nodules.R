###安装包
BiocManager::install("BioCor")
###
library("BioCor")
library("org.Hs.eg.db")
#library("reactome.db")
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(reshape2)
library(corrplot)
library(plyr)
library(igraph) 

###############################数据准备
###参考文献中用的GO数据库中BP分类下的通路，当然也可使使用其他通路
hsGO <- msigdbr(species = "Homo sapiens", category = "C5",subcategory='BP')
###差异分析结果
tumor_ec_markers0=read.xlsx('41467_2020_16164_MOESM8_ESM.xlsx',sheet = 3)

tumor_ec_markers0[1:5,]
Genes   log2FC       p-value    Bonferroni  PCT1  PCT2 
INSR 2.147045  7.910142e-57  1.328666e-52 0.585 0.126 
HSPG2 1.818224 4.165828e-120 6.997341e-116 0.864 0.321  
VWA1 1.716222  2.147301e-71  3.606822e-67 0.623 0.127       
PLVAP 1.695246  1.237309e-81  2.078309e-77 0.798 0.294  
IGHG3 1.694869  6.457962e-20  1.084744e-15 0.302 0.060    

###因为GO中很多通路是高度相似的，直接用会导致很多基因功能相似性过高，导致网络图中边的数量过多，所以先做富集得到重点通路
bp <- enrichGO(trio_up,keyType='SYMBOL', ont="BP", OrgDb = 'org.Hs.eg.db')
bp <- pairwise_termsim(bp)
bp2 <- clusterProfiler::simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)

###
GO_terms=bp2@result$ID#[bp2@result$pvalue<=0.05]
hsGO_tmp=hsGO[hsGO$gs_exact_source%in%GO_terms,]
hsGO_list=split(hsGO_tmp$gs_name,hsGO_tmp$gene_symbol)

###用Jaccard index系数计算基因之间的功能相似性
goSemSim <- mgeneSim(names(hsGO_list), info= hsGO_list,method = 'avg') ##avg , reciprocal
jaccard_index=D2J(goSemSim)

a=reshape2::melt(lower.tri(jaccard_index))
jaccard_index_df= reshape2::melt(jaccard_index)
jaccard_index_df=jaccard_index_df[a$value,]

#只保留值相似性较高的
jaccard_index_df0=jaccard_index_df[abs(jaccard_index_df$value)>0.15,]

############### 设置网络图的节点和边
nodes=trio_up
edges=jaccard_index_df0
nodes <- as.data.frame(nodes)
colnames(nodes) <- 'gene'
###设置节点的特征
#用pvalue控制节点的大小

#用FoldChange控制节点的颜色和深浅
nodes$cut=cut(nodes$log2FC,breaks = c(seq(min(nodes$log2FC)-0.1,-0.5,lengt=5),seq(0.5,max(nodes$log2FC)+0.1,lengt=5)))
# nodes$cut=cut(nodes$log2FC,10)
cols_nodes=colorRampPalette(colors = c("#143CFF","#F0F0F0","#FF0000"), interpolate ="linear")(9)
names(cols_nodes)=levels(nodes$cut)
nodes$color=cols_nodes[nodes$cut]

###设置边的特征
#edges$width=round(abs(edges$value)*10,0)
edges$width=c(2,2.5,8)[cut(edges$value,3)]
edges$color='#9C9C9C'

edges <- subset(edges,subset=Var1 %in% nodes$gene)
edges <- subset(edges,subset=Var2 %in% nodes$gene)
edges <-subset(edges,subset =value <3)
nodes <- subset(nodes,subset = gene %in% union(edges$Var1,edges$Var2))
###############开始画图
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 

V(net)$size <- 3
V(net)$frame.color <- '#4F4F4F'#nodes$frame_color
V(net)$shape='circle'
V(net)$label.color='black'
V(net)$label.size=0.1

E(net)$arrow.mode <- 0 #0:不需要箭头,1:back; 2:forth; 3:both
#E(net)$arrow.size=0.1 
#E(net)$arrow.width=0.1
E(net)$edge.color <- edges$color 
E(net)$width <- edges$width

plot(net,#vertex.label="",
     rescale = T,
     seed=1,
     layout=layout_nicely, 
     edge.curved=1
)

##############保存结果，以便用于Cytoscape画图
write.csv(nodes,'nodes.csv') 
write.csv(edges,'edges.csv')
