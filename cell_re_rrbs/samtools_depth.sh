data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs
output_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs/samtools_depth_result
folder=hepg2_bulk
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

cd ${data_path}/${folder}/sorted_bam
for bam in *.bam; do 
samtools depth -b ${bed} ${bam} > ${output_path}/${folder}/${bam%.bam}.depth
done

echo "done"