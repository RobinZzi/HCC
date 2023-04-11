cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/wes_hcc/hcc3/gatk_result/bqsr_result
  head -1 PT1.hg38_multianno.txt|sed 's/Otherinfo/Tumor_Sample_Barcode/' >header
cat header *maf > hcc3.maf