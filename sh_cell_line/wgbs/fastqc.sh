cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/sh_cell_line/wgbs
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs/qc
path=$1
files=$(ls$path)
for filename in $files
do
 cd $filename
 fastqc -t 12 -o $qcdir *.fastq.gz
 cd ..
done