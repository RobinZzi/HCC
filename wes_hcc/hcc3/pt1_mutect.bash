ref=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc3/wes_hcc3/bwa_result/hg38.exon.bed

    echo "start Mutect2 for PT1" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I BCPT1_bqsr.bam -tumor BCPT1 \
    -I BCNT_bqsr.bam -normal BCNT \
    -L $bed  \
    -O BCPT1_mutect2.vcf

    echo "end Mutect2 for PT1" `date`

