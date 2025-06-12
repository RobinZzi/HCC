bigmeth_colmeans_cmnsn <- bigmeth_colmeans


bigmeth_colmeans_cmnsn <- mutate(bigmeth_colmeans,state = case_when(origin %in% c("hcc3","hcc11","hcc28","hcc29")  ~ "CMN",
                                                                    origin %in% c("hcc7","hcc4") ~ "SN"))

bigmeth_colmeans_cmnsn$group <- paste(bigmeth_colmeans_cmnsn$state,bigmeth_colmeans_cmnsn$label,sep = "_")
bigmeth_colmeans_cmnsn$group <- factor(bigmeth_colmeans_cmnsn$group,levels = c("CMN_Normal Tissue","SN_Normal Tissue","CMN_Tumor Tissue","SN_Tumor Tissue"))



ggplot(bigmeth_colmeans_cmnsn,aes(x=group, y=mean_meth,fill=state))+
  geom_flat_violin(scale = "width",trim = F)+coord_flip()+geom_jitter(width = 0.1,size=0.5)+theme_bw()
