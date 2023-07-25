cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/oe
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/qc
path=$1
files=$(ls$path)
for filename in $files
do
 cd $filename
 fastqc -t 2 -o $qcdir *.fastq.gz
 cd ..
done

