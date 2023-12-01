HPC <- FindVariableFeatures(HPC, selection.method = "vst", nfeatures = 3000)
top3000_gene <- head(VariableFeatures(HPC), 3000)





hcc1_pt2 <- subset(HCC1_HPC, subset = orig.ident=="PT2")
hcc1_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt2, slot = "data"))
hcc1_pt2_top_mtx <- hcc1_pt2_mtx[_gene,]
hcc1_pt2_top_mtx <- na.omit(hcc1_pt2_top_mtx)
hcc1_pt2_cor <- as.data.frame(cor(hcc1_pt2_top_mtx,method = 'pearson'))


for(i in 1:nrow(hcc1_pt2_cor)){
  vector <- as.vector(hcc1_pt2_cor[i,1:(nrow(hcc1_pt2_cor)-i+1)])
  if(i == 1){
    hcc1_pt2_cor_nu <- vector
  }
  else{
    hcc1_pt2_cor_nu <- c(hcc1_pt2_cor_nu,vector)
  }
}



hcc1_pt4 <- subset(HCC1_HPC, subset = orig.ident=="PT4")
hcc1_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt4, slot = "data"))
hcc1_pt4_top_mtx <- hcc1_pt4_mtx[top3000_gene,]
hcc1_pt4_top_mtx <- na.omit(hcc1_pt4_top_mtx)
hcc1_pt4_cor <- as.data.frame(cor(hcc1_pt4_top_mtx,method = 'pearson'))
for(i in 1:nrow(hcc1_pt4_cor)){
  vector <- as.vector(hcc1_pt4_cor[i,1:(nrow(hcc1_pt4_cor)-i+1)])
  if(i == 1){
    hcc1_pt4_cor_nu <- vector
  }
  else{
    hcc1_pt4_cor_nu <- c(hcc1_pt4_cor_nu,vector)
  }
}

hcc1_pt5 <- subset(HCC1_HPC, subset = orig.ident=="PT5")
hcc1_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt5, slot = "data"))
hcc1_pt5_top_mtx <- hcc1_pt5_mtx[top3000_gene,]
hcc1_pt5_top_mtx <- na.omit(hcc1_pt5_top_mtx)
hcc1_pt5_cor <- as.data.frame(cor(hcc1_pt5_top_mtx,method = 'pearson'))
for(i in 1:nrow(hcc1_pt5_cor)){
  vector <- as.vector(hcc1_pt5_cor[i,1:(nrow(hcc1_pt5_cor)-i+1)])
  if(i == 1){
    hcc1_pt5_cor_nu <- vector
  }
  else{
    hcc1_pt5_cor_nu <- c(hcc1_pt5_cor_nu,vector)
  }
}


hcc1_primary_cor3 <- as.data.frame(c(hcc1_pt4_cor_nu,hcc1_pt2_cor_nu))
hcc1_primary_cor3 <-as.data.frame(t(hcc1_primary_cor3))
colnames(hcc1_primary_cor3) <- "cor"
hcc1_primary_cor3$type <- "primary"
hcc1_satellite_cor3 <- as.data.frame(unlist(unique(hcc1_pt5_cor_nu)))
colnames(hcc1_satellite_cor3) <- "cor"
hcc1_satellite_cor3$type <- "satellite"
hcc1_ps_info2 <- rbind(hcc1_primary_cor3,hcc1_satellite_cor3)



ggplot(data=hcc1_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC1")







hcc3_pt1 <- subset(HCC3_HPC, subset = orig.ident=="PT1")
hcc3_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt1, slot = "data"))
hcc3_pt1_top_mtx <- hcc3_pt1_mtx[top3000_gene,]
hcc3_pt1_top_mtx <- na.omit(hcc3_pt1_top_mtx)
hcc3_pt1_cor <- as.data.frame(cor(hcc3_pt1_top_mtx,method = 'pearson'))
hcc3_pt1_cor_nu <- unique(as.vector(as.matrix(hcc3_pt1_cor)))
for(i in 1:nrow(hcc3_pt1_cor)){
  vector <- as.vector(hcc3_pt1_cor[i,1:(nrow(hcc3_pt1_cor)-i+1)])
  if(i == 1){
    hcc3_pt1_cor_nu <- vector
  }
  else{
    hcc3_pt1_cor_nu <- c(hcc3_pt1_cor_nu,vector)
  }
}

hcc3_pt2 <- subset(HCC3_HPC, subset = orig.ident=="PT2")
hcc3_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt2, slot = "data"))
hcc3_pt2_top_mtx <- hcc3_pt2_mtx[top3000_gene,]
hcc3_pt2_top_mtx <- na.omit(hcc3_pt2_top_mtx)
hcc3_pt2_cor <- as.data.frame(cor(hcc3_pt2_top_mtx,method = 'pearson'))
hcc3_pt2_cor_nu <- unique(as.vector(as.matrix(hcc3_pt2_cor)))
for(i in 1:nrow(hcc3_pt2_cor)){
  vector <- as.vector(hcc3_pt2_cor[i,1:(nrow(hcc3_pt2_cor)-i+1)])
  if(i == 1){
    hcc3_pt2_cor_nu <- vector
  }
  else{
    hcc3_pt2_cor_nu <- c(hcc3_pt2_cor_nu,vector)
  }
}

hcc3_pt3 <- subset(HCC3_HPC, subset = orig.ident=="PT3")
hcc3_pt3_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt3, slot = "data"))
hcc3_pt3_top_mtx <- hcc3_pt3_mtx[top3000_gene,]
hcc3_pt3_top_mtx <- na.omit(hcc3_pt3_top_mtx)
hcc3_pt3_cor <- as.data.frame(cor(hcc3_pt3_top_mtx,method = 'pearson'))
hcc3_pt3_cor_nu <- unique(as.vector(as.matrix(hcc3_pt3_cor)))
for(i in 1:nrow(hcc3_pt3_cor)){
  vector <- as.vector(hcc3_pt3_cor[i,1:(nrow(hcc3_pt3_cor)-i+1)])
  if(i == 1){
    hcc3_pt3_cor_nu <- vector
  }
  else{
    hcc3_pt3_cor_nu <- c(hcc3_pt3_cor_nu,vector)
  }
}

hcc3_primary_cor3 <- as.data.frame(c(hcc3_pt1_cor_nu,hcc3_pt2_cor_nu))
hcc3_primary_cor3 <-as.data.frame(t(hcc3_primary_cor3))
colnames(hcc3_primary_cor3) <- "cor"
hcc3_primary_cor3$type <- "primary"
hcc3_satellite_cor3 <- as.data.frame(unlist(unique(hcc3_pt3_cor_nu)))
colnames(hcc3_satellite_cor3) <- "cor"
hcc3_satellite_cor3$type <- "satellite"
hcc3_ps_info2 <- rbind(hcc3_primary_cor3,hcc3_satellite_cor3)
ggplot(data=hcc3_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC3")


######hcc5-pt1/pt2/pt3/pt4####pt5#####
hcc5_pt1 <- subset(HCC5_HPC, subset = orig.ident=="PT1")
hcc5_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt1, slot = "data"))
hcc5_pt1_top_mtx <- hcc5_pt1_mtx[top3000_gene,]
hcc5_pt1_top_mtx <- na.omit(hcc5_pt1_top_mtx)
hcc5_pt1_cor <- as.data.frame(cor(hcc5_pt1_top_mtx,method = 'pearson'))
hcc5_pt1_cor_nu <- unique(as.vector(as.matrix(hcc5_pt1_cor)))
for(i in 1:nrow(hcc5_pt1_cor)){
  vector <- as.vector(hcc5_pt1_cor[i,1:(nrow(hcc5_pt1_cor)-i+1)])
  if(i == 1){
    hcc5_pt1_cor_nu <- vector
  }
  else{
    hcc5_pt1_cor_nu <- c(hcc5_pt1_cor_nu,vector)
  }
}




hcc5_pt2 <- subset(HCC5_HPC, subset = orig.ident=="PT2")
hcc5_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt2, slot = "data"))
hcc5_pt2_top_mtx <- hcc5_pt2_mtx[top3000_gene,]
hcc5_pt2_top_mtx <- na.omit(hcc5_pt2_top_mtx)
hcc5_pt2_cor <- as.data.frame(cor(hcc5_pt2_top_mtx,method = 'pearson'))
hcc5_pt2_cor_nu <- unique(as.vector(as.matrix(hcc5_pt2_cor)))
for(i in 1:nrow(hcc5_pt2_cor)){
  vector <- as.vector(hcc5_pt2_cor[i,1:(nrow(hcc5_pt2_cor)-i+1)])
  if(i == 1){
    hcc5_pt2_cor_nu <- vector
  }
  else{
    hcc5_pt2_cor_nu <- c(hcc5_pt2_cor_nu,vector)
  }
}






hcc5_pt3 <- subset(HCC5_HPC, subset = orig.ident=="PT3")
hcc5_pt3_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt3, slot = "data"))
hcc5_pt3_top_mtx <- hcc5_pt3_mtx[top3000_gene,]
hcc5_pt3_top_mtx <- na.omit(hcc5_pt3_top_mtx)
hcc5_pt3_cor <- as.data.frame(cor(hcc5_pt3_top_mtx,method = 'pearson'))
hcc5_pt3_cor_nu <- unique(as.vector(as.matrix(hcc5_pt3_cor)))
for(i in 1:nrow(hcc5_pt3_cor)){
  vector <- as.vector(hcc5_pt3_cor[i,1:(nrow(hcc5_pt3_cor)-i+1)])
  if(i == 1){
    hcc5_pt3_cor_nu <- vector
  }
  else{
    hcc5_pt3_cor_nu <- c(hcc5_pt3_cor_nu,vector)
  }
}




hcc5_pt4 <- subset(HCC5_HPC, subset = orig.ident=="PT4")
hcc5_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt4, slot = "data"))
hcc5_pt4_top_mtx <- hcc5_pt4_mtx[top3000_gene,]
hcc5_pt4_top_mtx <- na.omit(hcc5_pt4_top_mtx)
hcc5_pt4_cor <- as.data.frame(cor(hcc5_pt4_top_mtx,method = 'pearson'))
hcc5_pt4_cor_nu <- unique(as.vector(as.matrix(hcc5_pt4_cor)))
for(i in 1:nrow(hcc5_pt4_cor)){
  vector <- as.vector(hcc5_pt4_cor[i,1:(nrow(hcc5_pt4_cor)-i+1)])
  if(i == 1){
    hcc5_pt4_cor_nu <- vector
  }
  else{
    hcc5_pt4_cor_nu <- c(hcc5_pt4_cor_nu,vector)
  }
}



hcc5_pt5 <- subset(HCC5_HPC, subset = orig.ident=="PT5")
hcc5_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt5, slot = "data"))
hcc5_pt5_top_mtx <- hcc5_pt5_mtx[top3000_gene,]
hcc5_pt5_top_mtx <- na.omit(hcc5_pt5_top_mtx)
hcc5_pt5_cor <- as.data.frame(cor(hcc5_pt5_top_mtx,method = 'pearson'))
hcc5_pt5_cor_nu <- unique(as.vector(as.matrix(hcc5_pt5_cor)))
for(i in 1:nrow(hcc5_pt5_cor)){
  vector <- as.vector(hcc5_pt5_cor[i,1:(nrow(hcc5_pt5_cor)-i+1)])
  if(i == 1){
    hcc5_pt5_cor_nu <- vector
  }
  else{
    hcc5_pt5_cor_nu <- c(hcc5_pt5_cor_nu,vector)
  }
}



hcc5_primary_cor3 <- as.data.frame(c(hcc5_pt1_cor_nu,hcc5_pt2_cor_nu,hcc5_pt3_cor_nu,hcc5_pt4_cor_nu))
hcc5_primary_cor3 <-as.data.frame(t(hcc5_primary_cor3))
colnames(hcc5_primary_cor3) <- "cor"
hcc5_primary_cor3$type <- "primary"
hcc5_satellite_cor3 <- as.data.frame(unlist(unique(hcc5_pt5_cor_nu)))
colnames(hcc5_satellite_cor3) <- "cor"
hcc5_satellite_cor3$type <- "satellite"
hcc5_ps_info2 <- rbind(hcc5_primary_cor3,hcc5_satellite_cor3)
ggplot(data=hcc5_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC5")

######hcc6-pt1/pt2####pt4/pt5/pt6#####
hcc6_pt1 <- subset(HCC6_HPC, subset = orig.ident=="PT1")
hcc6_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt1, slot = "data"))
hcc6_pt1_top_mtx <- hcc6_pt1_mtx[top3000_gene,]
hcc6_pt1_top_mtx <- na.omit(hcc6_pt1_top_mtx)
hcc6_pt1_cor <- as.data.frame(cor(hcc6_pt1_top_mtx,method = 'pearson'))
hcc6_pt1_cor_nu <- unique(as.vector(as.matrix(hcc6_pt1_cor)))
for(i in 1:nrow(hcc6_pt1_cor)){
  vector <- as.vector(hcc6_pt1_cor[i,1:(nrow(hcc6_pt1_cor)-i+1)])
  if(i == 1){
    hcc6_pt1_cor_nu <- vector
  }
  else{
    hcc6_pt1_cor_nu <- c(hcc6_pt1_cor_nu,vector)
  }
}




hcc6_pt2 <- subset(HCC6_HPC, subset = orig.ident=="PT2")
hcc6_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt2, slot = "data"))
hcc6_pt2_top_mtx <- hcc6_pt2_mtx[top3000_gene,]
hcc6_pt2_top_mtx <- na.omit(hcc6_pt2_top_mtx)
hcc6_pt2_cor <- as.data.frame(cor(hcc6_pt2_top_mtx,method = 'pearson'))
hcc6_pt2_cor_nu <- unique(as.vector(as.matrix(hcc6_pt2_cor)))
for(i in 1:nrow(hcc6_pt2_cor)){
  vector <- as.vector(hcc6_pt2_cor[i,1:(nrow(hcc6_pt2_cor)-i+1)])
  if(i == 1){
    hcc6_pt2_cor_nu <- vector
  }
  else{
    hcc6_pt2_cor_nu <- c(hcc6_pt2_cor_nu,vector)
  }
}





hcc6_pt4 <- subset(HCC6_HPC, subset = orig.ident=="PT4")
hcc6_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt4, slot = "data"))
hcc6_pt4_top_mtx <- hcc6_pt4_mtx[top3000_gene,]
hcc6_pt4_top_mtx <- na.omit(hcc6_pt4_top_mtx)
hcc6_pt4_cor <- as.data.frame(cor(hcc6_pt4_top_mtx,method = 'pearson'))
hcc6_pt4_cor_nu <- unique(as.vector(as.matrix(hcc6_pt4_cor)))
for(i in 1:nrow(hcc6_pt4_cor)){
  vector <- as.vector(hcc6_pt4_cor[i,1:(nrow(hcc6_pt4_cor)-i+1)])
  if(i == 1){
    hcc6_pt4_cor_nu <- vector
  }
  else{
    hcc6_pt4_cor_nu <- c(hcc6_pt4_cor_nu,vector)
  }
}




hcc6_pt5 <- subset(HCC6_HPC, subset = orig.ident=="PT5")
hcc6_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt5, slot = "data"))
hcc6_pt5_top_mtx <- hcc6_pt5_mtx[top3000_gene,]
hcc6_pt5_top_mtx <- na.omit(hcc6_pt5_top_mtx)
hcc6_pt5_cor <- as.data.frame(cor(hcc6_pt5_top_mtx,method = 'pearson'))
hcc6_pt5_cor_nu <- unique(as.vector(as.matrix(hcc6_pt5_cor)))
for(i in 1:nrow(hcc6_pt5_cor)){
  vector <- as.vector(hcc6_pt5_cor[i,1:(nrow(hcc6_pt5_cor)-i+1)])
  if(i == 1){
    hcc6_pt5_cor_nu <- vector
  }
  else{
    hcc6_pt5_cor_nu <- c(hcc6_pt5_cor_nu,vector)
  }
}



hcc6_pt6 <- subset(HCC6_HPC, subset = orig.ident=="PT6")
hcc6_pt6_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt6, slot = "data"))
hcc6_pt6_top_mtx <- hcc6_pt6_mtx[top3000_gene,]
hcc6_pt6_top_mtx <- na.omit(hcc6_pt6_top_mtx)
hcc6_pt6_cor <- as.data.frame(cor(hcc6_pt6_top_mtx,method = 'pearson'))
hcc6_pt6_cor_nu <- unique(as.vector(as.matrix(hcc6_pt6_cor)))
for(i in 1:nrow(hcc6_pt6_cor)){
  vector <- as.vector(hcc6_pt6_cor[i,1:(nrow(hcc6_pt6_cor)-i+1)])
  if(i == 1){
    hcc6_pt6_cor_nu <- vector
  }
  else{
    hcc6_pt6_cor_nu <- c(hcc6_pt6_cor_nu,vector)
  }
}






hcc6_primary_cor3 <- as.data.frame(c(hcc6_pt1_cor_nu,hcc6_pt2_cor_nu))
hcc6_primary_cor3 <-as.data.frame(t(hcc6_primary_cor3))
colnames(hcc6_primary_cor3) <- "cor"
hcc6_primary_cor3$type <- "primary"
hcc6_satellite_cor3 <- as.data.frame(c(hcc6_pt4_cor_nu,hcc6_pt5_cor_nu,hcc6_pt6_cor_nu))
hcc6_satellite_cor3 <-as.data.frame(t(hcc6_satellite_cor3))
colnames(hcc6_satellite_cor3) <- "cor"
hcc6_satellite_cor3$type <- "satellite"
hcc6_ps_info2 <- rbind(hcc6_primary_cor3,hcc6_satellite_cor3)
ggplot(data=hcc6_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC6")





######hcc7-pt1/pt2/pt3/pt4/pt5###pt6##

hcc7_pt1 <- subset(HCC7_HPC, subset = orig.ident=="PT1")
hcc7_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt1, slot = "data"))
hcc7_pt1_top_mtx <- hcc7_pt1_mtx[top3000_gene,]
hcc7_pt1_top_mtx <- na.omit(hcc7_pt1_top_mtx)
hcc7_pt1_cor <- as.data.frame(cor(hcc7_pt1_top_mtx,method = 'pearson'))
hcc7_pt1_cor_nu <- unique(as.vector(as.matrix(hcc7_pt1_cor)))
for(i in 1:nrow(hcc7_pt1_cor)){
  vector <- as.vector(hcc7_pt1_cor[i,1:(nrow(hcc7_pt1_cor)-i+1)])
  if(i == 1){
    hcc7_pt1_cor_nu <- vector
  }
  else{
    hcc7_pt1_cor_nu <- c(hcc7_pt1_cor_nu,vector)
  }
}



hcc7_pt2 <- subset(HCC7_HPC, subset = orig.ident=="PT2")
hcc7_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt2, slot = "data"))
hcc7_pt2_top_mtx <- hcc7_pt2_mtx[top3000_gene,]
hcc7_pt2_top_mtx <- na.omit(hcc7_pt2_top_mtx)
hcc7_pt2_cor <- as.data.frame(cor(hcc7_pt2_top_mtx,method = 'pearson'))
hcc7_pt2_cor_nu <- unique(as.vector(as.matrix(hcc7_pt2_cor)))
for(i in 1:nrow(hcc7_pt2_cor)){
  vector <- as.vector(hcc7_pt2_cor[i,1:(nrow(hcc7_pt2_cor)-i+1)])
  if(i == 1){
    hcc7_pt2_cor_nu <- vector
  }
  else{
    hcc7_pt2_cor_nu <- c(hcc7_pt2_cor_nu,vector)
  }
}



hcc7_pt3 <- subset(HCC7_HPC, subset = orig.ident=="PT3")
hcc7_pt3_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt3, slot = "data"))
hcc7_pt3_top_mtx <- hcc7_pt3_mtx[top3000_gene,]
hcc7_pt3_top_mtx <- na.omit(hcc7_pt3_top_mtx)
hcc7_pt3_cor <- as.data.frame(cor(hcc7_pt3_top_mtx,method = 'pearson'))
hcc7_pt3_cor_nu <- unique(as.vector(as.matrix(hcc7_pt3_cor)))
for(i in 1:nrow(hcc7_pt3_cor)){
  vector <- as.vector(hcc7_pt3_cor[i,1:(nrow(hcc7_pt3_cor)-i+1)])
  if(i == 1){
    hcc7_pt3_cor_nu <- vector
  }
  else{
    hcc7_pt3_cor_nu <- c(hcc7_pt3_cor_nu,vector)
  }
}




hcc7_pt4 <- subset(HCC7_HPC, subset = orig.ident=="PT4")
hcc7_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt4, slot = "data"))
hcc7_pt4_top_mtx <- hcc7_pt4_mtx[top3000_gene,]
hcc7_pt4_top_mtx <- na.omit(hcc7_pt4_top_mtx)
hcc7_pt4_cor <- as.data.frame(cor(hcc7_pt4_top_mtx,method = 'pearson'))
hcc7_pt4_cor_nu <- unique(as.vector(as.matrix(hcc7_pt4_cor)))
for(i in 1:nrow(hcc7_pt4_cor)){
  vector <- as.vector(hcc7_pt4_cor[i,1:(nrow(hcc7_pt4_cor)-i+1)])
  if(i == 1){
    hcc7_pt4_cor_nu <- vector
  }
  else{
    hcc7_pt4_cor_nu <- c(hcc7_pt4_cor_nu,vector)
  }
}




hcc7_pt5 <- subset(HCC7_HPC, subset = orig.ident=="PT5")
hcc7_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt5, slot = "data"))
hcc7_pt5_top_mtx <- hcc7_pt5_mtx[top3000_gene,]
hcc7_pt5_top_mtx <- na.omit(hcc7_pt5_top_mtx)
hcc7_pt5_cor <- as.data.frame(cor(hcc7_pt5_top_mtx,method = 'pearson'))
hcc7_pt5_cor_nu <- unique(as.vector(as.matrix(hcc7_pt5_cor)))
for(i in 1:nrow(hcc7_pt5_cor)){
  vector <- as.vector(hcc7_pt5_cor[i,1:(nrow(hcc7_pt5_cor)-i+1)])
  if(i == 1){
    hcc7_pt5_cor_nu <- vector
  }
  else{
    hcc7_pt5_cor_nu <- c(hcc7_pt5_cor_nu,vector)
  }
}






hcc7_pt6 <- subset(HCC7_HPC, subset = orig.ident=="PT6")
hcc7_pt6_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt6, slot = "data"))
hcc7_pt6_top_mtx <- hcc7_pt6_mtx[top3000_gene,]
hcc7_pt6_top_mtx <- na.omit(hcc7_pt6_top_mtx)
hcc7_pt6_cor <- as.data.frame(cor(hcc7_pt6_top_mtx,method = 'pearson'))
hcc7_pt6_cor_nu <- unique(as.vector(as.matrix(hcc7_pt6_cor)))
for(i in 1:nrow(hcc7_pt6_cor)){
  vector <- as.vector(hcc7_pt6_cor[i,1:(nrow(hcc7_pt6_cor)-i+1)])
  if(i == 1){
    hcc7_pt6_cor_nu <- vector
  }
  else{
    hcc7_pt6_cor_nu <- c(hcc7_pt6_cor_nu,vector)
  }
}






hcc7_primary_cor3 <- as.data.frame(c(hcc7_pt1_cor_nu,hcc7_pt2_cor_nu,hcc7_pt3_cor_nu,hcc7_pt4_cor_nu,hcc7_pt5_cor_nu))
hcc7_primary_cor3 <-as.data.frame(t(hcc7_primary_cor3))
colnames(hcc7_primary_cor3) <- "cor"
hcc7_primary_cor3$type <- "primary"
hcc7_satellite_cor3 <- as.data.frame(unlist(unique(hcc7_pt6_cor_nu)))
colnames(hcc7_satellite_cor3) <- "cor"
hcc7_satellite_cor3$type <- "satellite"
hcc7_ps_info2 <- rbind(hcc7_primary_cor3,hcc7_satellite_cor3)
ggplot(data=hcc7_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC7")

hcc1_ps_info2$sample <- "hcc1"
hcc3_ps_info2$sample <- "hcc3"
hcc5_ps_info2$sample <- "hcc5"
hcc6_ps_info2$sample <- "hcc6"
hcc7_ps_info2$sample <- "hcc7"

sta_prim_sum <- rbind(hcc1_ps_info2,hcc3_ps_info2,hcc5_ps_info2,
                      hcc6_ps_info2,hcc7_ps_info2)
ggplot(sta_prim_sum,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("primary","satellite")),
                     method = "t.test",label = "p.signif",
                     label.y =1)+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())
ggplot(sta_prim_sum,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~sample)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("primary","satellite")),
                     method = "t.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())
