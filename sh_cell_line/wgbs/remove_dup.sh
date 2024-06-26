dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240611/bismark_bam
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240611/remove_dup

cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
deduplicate_bismark --bam ${headname}.bam --output_dir $dir2
echo $headname
done
