cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex3/trim_result
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex3/aftertrim_qc

fastqc -t 8 -o $qcdir *.fq.gz
