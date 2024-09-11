data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs
ref_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs/cnvkit_result/ref
geno_path=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
output_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs/cnvkit_result/result
GENOME=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed
folder=hepg2_trioseq_2

cd ${data_path}/${folder}/sorted_bam

for bam in *.bam; do 
cnvkit.py batch ${bam} \
 -r ${ref_path}/my_reference.cnn\
 -p 10 \
 --drop-low-coverage --scatter --diagram --method amplicon \
 --output-dir ${output_path}/${folder}/ 
done


echo '##############################'
echo '############done##############'
echo '##############################'