ref=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/wes_hcc/hcc3/gatk_result/bqsr_result
for sample in {PT1,PT2,PT3,PT4}
do
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  FilterMutectCalls -R $ref\
    -V ${sample}_mutect2.vcf \
    -O ${sample}_somatic.vcf 
    echo "end Mutect2 for ${sample}" `date`
done  