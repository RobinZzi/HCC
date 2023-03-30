ref=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed

cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc29/gatk_result
    echo "start Mutect2 for PT3" `date`
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I PT3_bqsr.bam -tumor PT3 \
    -I NT_bqsr.bam -normal NT \
    -L $bed  \
    -O PT3_mutect2.vcf

    echo "end Mutect2 for PT3" `date`
    
    
    
