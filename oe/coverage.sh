dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/remove_dup
outP=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/coverage
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38


cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
bismark_methylation_extractor -p --no_overlap --parallel 30 --buffer_size 90G --counts --output ${outP} --genome_folder ${index} ${headname}.bam
done
