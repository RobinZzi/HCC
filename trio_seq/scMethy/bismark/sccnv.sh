bam_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/sorted_bam
ref_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/cnvkit_result/scMerge/ref
geno_path=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
output_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/cnvkit_result/scMerge
GENOME=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

#cd ${bam_path}/nt

#for bam in nt*.bam; do  
#    cnvkit.py coverage $bam ${bed} -o ${ref_path}/${bam%.bam}targetcoverage.cnn  
#done

echo '##############################'
echo '##########nt prepared#########'
echo '##############################'

#cd ${ref_path}
#cnvkit.py reference *.cnn -f $GENOME -o my_reference.cnn

echo '##############################'
echo '#######ref_sum prepared#######'
echo '##############################'


cd ${bam_path}
for pt in pt*; do
cd ${pt}
for bam in pt*.bam; do 
cd ${bam_path}/${pt}
cnvkit.py batch pt*_sorted.bam \
 -r ${ref_path}/my_reference.cnn\
 -p 10 \
 --drop-low-coverage --scatter --diagram --method amplicon \
 --output-dir ${output_path}/${pt}/ 
done
done

echo '##############################'
echo '############done##############'
echo '##############################'