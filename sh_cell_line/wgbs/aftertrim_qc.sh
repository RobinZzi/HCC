cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs/trim_result
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs/aftertrim_qc

fastqc -t 8 -o $qcdir *.fq.gz