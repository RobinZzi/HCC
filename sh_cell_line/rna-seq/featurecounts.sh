dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/featurecounts_result
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/samtools_result
ref=/storage/zhangyanxiaoLab/zhangliwen/genome/gencode.v40.annotation.gtf.gz

    #featureCounts -T 16 -a $ref -o ${dir}/matrix.txt *.bam
    
  #  ls *.bam |while read id
  #  do
  #  headname=${id%.bam}
  #  echo ${headname}
  #  featureCounts -T 16 -p -a $ref -o ${dir}/${headname}.count ${headname}.bam
  #  done

#featureCounts -T 16 -a $ref -t gene -f -g gene_id --extraAttributes gene_name -M -o ${dir}/matrix.txt *.bam
featureCounts -a genes.saf -o ${dir}/matrix.txt *.bam -F SAF -T 8