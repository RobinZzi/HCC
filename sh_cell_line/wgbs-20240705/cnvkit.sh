GENOME=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240705/bismark_bam
cnvkit.py batch sorted_CR2.bam \
--normal  /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_bam/sorted_n.bam \
--targets  ${bed} \
--fasta $GENOME  \
--drop-low-coverage --scatter --diagram --method amplicon \
--output-reference my_reference.cnn --output-dir /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240705/cnv/