cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/liver_cancer_cellline/chip_seq

#bigWigToBedGraph GSM2360939_1_ctrl_untr_Huh7_H3K27ac.bw hg19_Huh7_H3K27ac.bedGraph
#bigWigToBedGraph GSM2360948_1_ctrl_untr_Huh7_H3K9ac.bw hg19_Huh7_H3K9ac.bedGraph
#bigWigToBedGraph GSM2360945_1_ctrl_untr_Huh7_H3K4me1.bw hg19_Huh7_H3K4me1.bedGraph
bigWigToBedGraph GSM2360942_4_ctrl_untr_Huh7_H3K36ac.bw hg19_Huh7_H3K36ac.bedGraph
bigWigToBedGraph GSM2360951_4_ctrl_untr_Huh7_H4K5ac.bw hg19_Huh7_H4K5ac.bedGraph
 
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H3K27ac.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H3K27ac.bedGraph unmapped.bed  
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H3K9ac.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H3K9ac.bedGraph unmapped.bed
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H3K4me1.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H3K4me1.bedGraph unmapped.bed
#/storage/zhangyanxiaoLab/zhangliwen/src/liftOver GSE89212_enhancers.bed /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_huh7_enhancer.bed unmapped.bed


/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H4K5ac.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H4K5ac.bedGraph unmapped.bed  
/storage/zhangyanxiaoLab/zhangliwen/src/liftOver hg19_Huh7_H3K36ac.bedGraph /storage/zhangyanxiaoLab/zhangliwen/genome/hg19ToHg38.over.chain.gz hg38_Huh7_H3K36ac.bedGraph unmapped.bed  

