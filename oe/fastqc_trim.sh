cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex2/trim_result
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex2/aftertrim_qc

fastqc -t 8 -o $qcdir *.fq.gz
