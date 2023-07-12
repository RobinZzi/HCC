GENOME=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
bed=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/hg38.exon.bed
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc29/gatk_result/
cnvkit.py batch PT*bqsr.bam \
--normal  NT_bqsr.bam \
--targets  ${bed} \
--fasta $GENOME  \
--drop-low-coverage --scatter --diagram --method amplicon \
--output-reference my_reference.cnn --output-dir /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc29/gatk_result/cnvkit_result
