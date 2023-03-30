cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc28/gatk_result/
head -1 PT1.hg38_multianno.txt|sed 's/Otherinfo/Tumor_Sample_Barcode/' >header
cat header *maf > hcc28.maf
