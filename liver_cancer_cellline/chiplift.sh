cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/liver_cancer_cellline/chip_seq

#bigWigToBedGraph GSM2360939_1_ctrl_untr_Huh7_H3K27ac.bw hg19_Huh7_H3K27ac.bedGraph
#bigWigToBedGraph GSM2360948_1_ctrl_untr_Huh7_H3K9ac.bw hg19_Huh7_H3K9ac.bedGraph
#bigWigToBedGraph GSM2360945_1_ctrl_untr_Huh7_H3K4me1.bw hg19_Huh7_H3K4me1.bedGraph
bigWigToBedGraph GSM1249889_H3K27ac-HEK293T.bigWig hg18_293T_H3K27ac.bedGraph
bigWigToBedGraph GSM1249885_H3K4me3-HEK293T.bigWig hg18_293T_H3K4me3.bedGraph
 
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H3K27ac.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H3K27ac.bedGraph unmapped.bed  
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H3K9ac.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H3K9ac.bedGraph unmapped.bed
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H3K4me1.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H3K4me1.bedGraph unmapped.bed
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver GSE89212_enhancers.bed /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_huh7_enhancer.bed unmapped.bed


/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg18_293T_H3K27ac.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg18ToHg38.over.chain.gz hg38_293T_H3K27ac.bedGraph unmapped.bed  
/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg18_293T_H3K4me3.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg18ToHg38.over.chain.gz hg38_293T_H3K4me3.bedGraph unmapped.bed  

bigWigToBedGraph GSM7892792_HEK293T_H3K4me3_WT4_treat_pileup.bw hg38_293T_H3K4me3_2.bedGraph
bigWigToBedGraph GSM7892793_HEK293T_H3K9me3_WT9_treat_pileup.bw hg38_293T_H3K9me3.bedGraph