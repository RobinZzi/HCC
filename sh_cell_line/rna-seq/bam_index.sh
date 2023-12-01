dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/samtools_result
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/samtools_result
ls *.bam |while read id
    do
    headname=${id%.bam}
    echo ${headname}
    samtools index ${headname}.bam
    samtools flagstat ${headname}.bam
    done