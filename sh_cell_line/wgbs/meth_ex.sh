dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/remove_dup
outP=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/meth_ex
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38


cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
bismark_methylation_extractor -p --no_overlap --parallel 30 --bedGraph --buffer_size 90G --cytosine_report --output ${outP} --genome_folder ${index} ${headname}.bam
done

