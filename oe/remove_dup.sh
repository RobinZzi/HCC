
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/bismark_bam
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/remove_dup

cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
deduplicate_bismark --bam ${headname}.bam --o $dir
echo $headname
mv bismark_bam.deduplicated.bam $dir2/${headname}remove_dup.bam
cd $dir
done

