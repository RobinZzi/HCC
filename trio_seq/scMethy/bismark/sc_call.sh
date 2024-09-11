bam_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/sorted_bam
ref_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/cnvkit_result/scMerge/ref
output_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/cnvkit_result/scMerge
pt=pt1
for bam in pt*.bam; do 
cd ${bam_path}/${pt}
cnvkit.py batch pt*_sorted.bam \
 -r ${ref_path}/my_reference.cnn\
 -p 15 \
 --drop-low-coverage --scatter --diagram --method amplicon \
 --output-dir ${output_path}/${pt}/ 
 done
 