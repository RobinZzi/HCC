GENOME=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/merged_bam/sorted_bam
cnvkit.py batch pt*_sorted.bam \
--normal  nt_sorted.bam \
--targets  ${bed} \
--fasta $GENOME  \
--drop-low-coverage --scatter --diagram --method amplicon \
--output-reference my_reference.cnn\ 
--output-dir /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/cnvkit_result/merged_bam


