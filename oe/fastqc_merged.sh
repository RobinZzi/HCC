cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/cated_fastq
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/qc

fastqc -t 8 -o $qcdir *.fq.gz