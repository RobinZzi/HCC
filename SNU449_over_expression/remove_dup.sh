
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/over_expression/remove_dup
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/over_expression/bismark_bam

ls *.bam |while read id
do
headname=${id%.bam}
deduplicate_bismark --bam ${headname}.bam --o $dir
echo $headname
mv remove_dup.deduplicated.bam $dir
cd $dir
pwd
mv remove_dup.deduplicated.bam ${headname}_remove_dup.bam
cd $dir2
pwd
done
