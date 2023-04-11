ref=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/wes_hcc/hcc3/gatk_result/bqsr_result
    echo "start Mutect2 for PT2" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I BCPT1_bqsr.bam -tumor BCPT1 \
    -I BCNT_bqsr.bam -normal BCNT \
    -L $bed  \
    -O PT1_mutect2.vcf

    echo "end Mutect2 for PT1" `date`
    echo "start Mutect2 for PT2" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I BCPT2_bqsr.bam -tumor BCPT2 \
    -I BCNT_bqsr.bam -normal BCNT \
    -L $bed  \
    -O PT2_mutect2.vcf

    echo "end Mutect2 for PT2" `date`
    echo "start Mutect2 for PT3" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I BCPT3_bqsr.bam -tumor BCPT3 \
    -I BCNT_bqsr.bam -normal BCNT \
    -L $bed  \
    -O PT3_mutect2.vcf

    echo "end Mutect2 for PT3" `date`
    echo "start Mutect2 for PT4" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I BCPT4_bqsr.bam -tumor BCPT3 \
    -I BCNT_bqsr.bam -normal BCNT \
    -L $bed  \
    -O PT4_mutect2.vcf

    echo "end Mutect2 for PT4" `date`