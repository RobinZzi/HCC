data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs
ref_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs/cnvkit_result/ref
geno_path=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
output_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs/cnvkit_result/result
GENOME=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

cd ${data_path}/liver_bulk/sorted_bam
for bam in *.bam; do  
    cnvkit.py coverage $bam ${bed} -o ${ref_path}/${bam%.bam}targetcoverage.cnn  
done

echo '##############################'
echo '##########nt prepared#########'
echo '##############################'

cd ${ref_path}
cnvkit.py reference *.cnn -f $GENOME -o my_reference.cnn

echo '##############################'
echo '#######ref_sum prepared#######'
echo '##############################'