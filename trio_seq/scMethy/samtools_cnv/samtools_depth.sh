data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/sorted_bam
output_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/samtools_depth_result
folder=nt
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

cd ${data_path}/${folder}
for bam in *.bam; do 
samtools depth -b ${bed} ${bam} > ${output_path}/${folder}/${bam%.bam}.depth
done

echo '##############################'
echo "Current date and time: $(date)"
echo '############done##############'
echo '##############################'