cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/oe/ex3
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/ex3/qc
path=$1
files=$(ls$path)
for filename in $files
do
 cd $filename
 fastqc -t 12 -o $qcdir *.fastq.gz
 cd ..
done

