ls *bam | while read id
do
samtools stats -@ 16 --reference //storage/zhangyanxiaoLab/share/fasta/hg38.fa ${id} > /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc3/wes_hcc3/bwa_result/stat/$(basename ${id} .bam).stat
done
