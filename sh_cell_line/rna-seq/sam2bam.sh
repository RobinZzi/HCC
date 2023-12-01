dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/samtools_result
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/hisat2_result
ls *.SAM |while read id
    do
    headname=${id%.SAM}
    samtools view  -@ 4 -bS  ${headname}.SAM  >  ${headname}_unsorted.bam
    samtools sort ${headname}_unsorted.bam -o ${dir}/${headname}_sorted.bam
    done
    
