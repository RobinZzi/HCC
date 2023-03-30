ref=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc28/gatk_result
    echo "start Mutect2 for PT1" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I PT1_bqsr.bam -tumor PT1 \
    -I NT_bqsr.bam -normal NT \
    -L $bed  \
    -O PT1_mutect2.vcf

    echo "end Mutect2 for PT1" `date`
