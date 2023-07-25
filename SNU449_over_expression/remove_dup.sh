dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/over_expression/remove_dup
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/over_expression/bismark_bam

cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
deduplicate_bismark --bam ${headname}.bam --o $dir
mv bismark_bam.deduplicated.bam ${headname}_remove_dup.bam
cd $dir2
pwd
done
