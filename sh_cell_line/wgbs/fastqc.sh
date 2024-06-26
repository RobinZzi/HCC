cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240614/fastq
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240614/qc
path=$1
echo $path
files=$(ls$path)
for filename in $files
do
 cd $filename
 fastqc -t 12 -o $qcdir *.fastq.gz
 cd ..
done