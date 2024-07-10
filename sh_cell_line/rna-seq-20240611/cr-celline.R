rm(list = ls())
  library(FactoMineR)
library(data.table)
library(dplyr)
library(edgeR)

setwd("~/projects/hcc/analysis/liver_cancer_cellline")
setwd("~/projects/hcc/analysis/293T_rna")
count2 <- fread("293T_te_counts.txt")
rpkm2 <- fread("combined-chrM.rpkm")


setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-20240611")

count <- fread("combined-chrM.counts")
rpkm <- fread("combined-chrM.rpkm")

setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-20240611/diff_analysis")
save.image("sh-celline.Rdata")








sum_rpkm <- as.data.frame(cbind(dplyr::select(rpkm,1,7,8,9,10),rpkm2[,7:10]))
rownames(sum_rpkm) <- sum_rpkm[,1]
sum_rpkm <- sum_rpkm[,-1]
colnames(sum_rpkm) <- c("v-2","v-1","cr-2","cr-1","wt-1","wt-2","wt-3","wt-4")


sum_count <- as.data.frame(dplyr::select(count,1,7,8,9,10))
rownames(sum_count) <- sum_count[,1]
sum_count <- sum_count[,-1]
colnames(sum_count) <- c("v-2","v-1","cr-2","cr-1")



genelist <- c("SNHG6","MET","GADD45A","ONECUT2")

submtx <- sum_rpkm[genelist,]
submtx <- t(submtx)


genelist2 <- c("DNMT1","DNMT3B","DNMT3A","TET1","TET2","TDG")
submtx2 <- sum_rpkm[genelist2,]
submtx2 <- t(submtx2)



group = c("sh","sh","v","v")
design = model.matrix(~0+group)


y = DGEList(counts=sum_count)
y<-calcNormFactors(y)
y<-estimateCommonDisp(y, rowsum.filter=5)
y<-estimateGLMTagwiseDisp(y,design) #
fit_tag<-glmFit(y,design)
lrt.tagwise<-glmLRT(fit_tag,contrast=c(1,-1))
pvals_tag <- lrt.tagwise$table$PValue
FDR_tag<- p.adjust(pvals_tag, method="BH")
out1 = cbind(GeneID=rownames(y),lrt.tagwise$table,FDR_tag)
out1 = out1[order(out1$FDR_tag),]
table(out1$FDR_tag<0.05)
sig1 = out1[which(out1$FDR_tag<0.05),]
setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-20240611/diff_analysis")
write.csv(sig1,"cr-v-diff.csv")


cr_v_diff_all <- out1 %>%
  mutate(expression = case_when(logFC <= -1 & FDR_tag < 0.05 ~ "Up-regulated", # 上调
                                logFC >= 1 & FDR_tag < 0.05 ~ "Down-regulated", # 下调
                                TRUE ~ "Unchanged")) # 不变



ggplot(cr_v_diff_all, aes(logFC, -log10(FDR_tag))) +
  geom_point(size = 0.4, aes(color = expression)) + # 根据expression水平进行着色
  xlab(expression("log"[2]*" fold change")) + # 修饰x轴题目
  ylab(expression("-log"[10]*" FDR")) + # 修饰y轴题目
  scale_x_continuous(limits = c(-15, 15)) +
  scale_color_manual(values = c("steelblue", "grey", "red"))+ # 添加三种颜色
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  labs(title="cr-GADD45A vs Vector diff gene")





setwd("~/projects/hcc/analysis/sh_cell_line/rna-seq-20240611/diff_analysis")
cr_v_diff <- fread("cr-v-diff.csv")
cr_v_up <- subset(cr_v_diff,subset = logFC<0)
cr_v_down <- subset(cr_v_diff,subset = logFC>0)

cr_v_up_go <- enrichGO(gene  = cr_v_up$GeneID,
                       OrgDb      = org.Hs.eg.db,
                       keyType    = 'SYMBOL',
                       ont        = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
cr_v_up_go <- as.data.frame(cr_v_up_go@result)
cr_v_up_go [,"logp"] <- -log10(cr_v_up_go$pvalue)
cr_v_up_go[,11:12] <- as.numeric(str_split_fixed(cr_v_up_go$GeneRatio,"/",2))
cr_v_up_go$GeneRatio <- cr_v_up_go[,11]/cr_v_up_go[,12]
ggplot(data = cr_v_up_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="Up-Regulated Pathways in cr-GADD45A 293T")



cr_v_down_go <- enrichGO(gene  = cr_v_down$GeneID,
                         OrgDb      = org.Hs.eg.db,
                         keyType    = 'SYMBOL',
                         ont        = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
cr_v_down_go <- as.data.frame(cr_v_down_go@result)
cr_v_down_go [,"logp"] <- -log10(cr_v_down_go$pvalue)
cr_v_down_go[,11:12] <- as.numeric(str_split_fixed(cr_v_down_go$GeneRatio,"/",2))
cr_v_down_go$GeneRatio <- cr_v_down_go[,11]/cr_v_down_go[,12]
ggplot(data = cr_v_down_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(Count),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="down-Regulated Pathways in cr-GADD45A 293T")






rpkm_t <- t(sum_rpkm)
rpkm_group <- data.frame(Sample = rownames(rpkm_t), Group = c("sh","sh","v","v","wt","wt","wt","wt"))
rpkm_pca <- PCA(rpkm_t,scale.unit = TRUE,graph = FALSE)
rpkm_pca_sample <- data.frame(rpkm_pca$ind$coord[ ,1:2])
rpkm_pca_sample$Sample=row.names(rpkm_pca_sample)
rpkm_pca_eig1 <- round(rpkm_pca$eig[1,2], 2)
rpkm_pca_eig2 <- round(rpkm_pca$eig[2,2],2 )
rpkm_pca_sample <- merge(rpkm_pca_sample,rpkm_group,by="Sample")
head(rpkm_pca_sample)
ggplot(data = rpkm_pca_sample, aes(x = Dim.1, y = Dim.2))+
  geom_point(aes(color = Group), size = 5) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),legend.key = element_rect(fill = 'transparent')) + 
  labs(x = paste('PCA1:', rpkm_pca_eig1, '%'), y = paste('PCA2:', rpkm_pca_eig2, '%'), color = '') + 
  scale_fill_manual(values = c('orange', 'purple'))+
  labs(title="SNHG6_pca")


cr_v_down_list <- cr_v_down$GeneID
cr_v_up_list <- cr_v_up$GeneID


promoter_meth <- fread("dmr_sig_promoter.bedGraph")

promoter_demeth <- subset(promoter_meth,subset= V15>0)
promoter_admeth <- subset(promoter_meth,subset= V15<0)

promoter_demeth_list <- unique(promoter_demeth$V5)
promoter_admeth_list <- unique(promoter_admeth$V5)


demethlist <- list(promoter_demeth_list,cr_v_up_list)

admethlist <- list(promoter_admeth_list,cr_v_up_list)



venn.diagram(demethlist, filename = 'demeth.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, ,category.names = c("promoter_demeth" , "rna_up"),
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex =0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB'),cex = 2.5 )


venn.diagram(admethlist, filename = 'admeth.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, ,category.names = c("promoter_admeth" , "rna_down"),
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex =0, cat.fontfamily = 'serif',vjust=5,
             col = c('#4D157D', '#84C7DB'),cex = 2.5 )



promoter_demeth_go <- enrichGO(gene  = promoter_demeth_list,
                               OrgDb      = org.Hs.eg.db,
                               keyType    = 'SYMBOL',
                               ont        = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05)
promoter_demeth_go <- as.data.frame(promoter_demeth_go@result)
promoter_demeth_go [,"logp"] <- -log10(promoter_demeth_go$pvalue)
promoter_demeth_go [,"LogFDR"] <- -log10(promoter_demeth_go$p.adjust)
promoter_demeth_go_enrichment_fold <- apply(promoter_demeth_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
promoter_demeth_go$EF <- promoter_demeth_go_enrichment_fold
ggplot(data = promoter_demeth_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogP", 
       y= " ",
       title="demethylated genes in 293T-GA45KD Pathways")






promoter_admeth_go <- enrichGO(gene  = promoter_admeth_list,
                               OrgDb      = org.Hs.eg.db,
                               keyType    = 'SYMBOL',
                               ont        = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05)
promoter_admeth_go <- as.data.frame(promoter_admeth_go@result)
promoter_admeth_go [,"logp"] <- -log10(promoter_admeth_go$pvalue)
promoter_admeth_go [,"LogFDR"] <- -log10(promoter_admeth_go$p.adjust)
promoter_admeth_go_enrichment_fold <- apply(promoter_admeth_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
promoter_admeth_go$EF <- promoter_admeth_go_enrichment_fold
ggplot(data = promoter_admeth_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogP", 
       y= " ",
       title="admethylated genes in 293T-GA45KD Pathways")


write.table(promoter_demeth_go,"promoter_demeth_go.csv",row.names = T,col.names = T,sep = "\t",quote=F)
write.table(promoter_admeth_go,"promoter_admeth_go.csv",row.names = T,col.names = T,sep = "\t",quote=F)
write.table(cr_v_down_go,"down_kd_go.csv",row.names = T,col.names = T,sep = "\t",quote=F)
write.table(cr_v_up_go,"up_kd_go",row.names = T,col.names = T,sep = "\t",quote=F)


promoter_admeth_go_sig <- subset(promoter_admeth_go,subset = promoter_admeth_go$pvalue <0.05)
promoter_demeth_go_sig <- subset(promoter_demeth_go,subset = promoter_demeth_go$pvalue <0.05)


cr_v_down_go_sig <- subset(cr_v_down_go,subset=cr_v_down_go$pvalue<0.05)
cr_v_up_go_sig <- subset(cr_v_up_go,subset=cr_v_up_go$pvalue<0.05)

cr_v_down_go_sig_2 <- subset(cr_v_down_go,subset=cr_v_down_go$p.adjust<0.05)
cr_v_up_go_sig_2 <- subset(cr_v_up_go,subset=cr_v_up_go$p.adjust<0.05)


intersect(promoter_admeth_go_sig$Description,cr_v_down_go_sig$Description)

intersect(promoter_demeth_go_sig$Description,cr_v_up_go_sig$Description)



intersect(promoter_admeth_go_sig$Description,cr_v_down_go_sig_2$Description)

intersect(promoter_demeth_go_sig$Description,cr_v_up_go_sig_2$Description)



demeth_pathway_list <- list(promoter_demeth_go_sig$Description,cr_v_up_go_sig$Description)

admeth_pathway_list <- list(promoter_admeth_go_sig$Description,cr_v_down_go_sig$Description)



venn.diagram(demeth_pathway_list, filename = 'demeth_pathway.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, ,category.names = c("promoter_demeth" , "rna_up"),
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex =0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB'),cex = 2.5 )


venn.diagram(admeth_pathway_list, filename = 'admeth_pathway.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, ,category.names = c("promoter_admeth" , "rna_down"),
             fill = c('#4D157D', '#84C7DB'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB'), cat.cex =0, cat.fontfamily = 'serif',vjust=5,
             col = c('#4D157D', '#84C7DB'),cex = 2.5 )

