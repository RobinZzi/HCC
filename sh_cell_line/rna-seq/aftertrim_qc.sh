cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/trim_result
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/aftertrim_qc

fastqc -t 8 -o $qcdir *.fq.gz