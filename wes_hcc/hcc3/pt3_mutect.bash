ref=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc3/wes_hcc3/bwa_result/hg38.exon.bed

    echo "start Mutect2 for PT3" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I BCPT3_bqsr.bam -tumor BCPT3 \
    -I BCNT_bqsr.bam -normal BCNT \
    -L $bed  \
    -O BCPT3_mutect2.vcf

    echo "end Mutect2 for PT3" `date`

