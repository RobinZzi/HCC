cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/oe
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/trim_result
path=$1
files=$(ls$path)
for filename in $files
do
 cd $filename
 ls *fastq.gz|while read id
 headname=${id%.R*}
 fq1=${headname}.R1.fastq.gz
 fq2=${headname}.R2.fastq.gz
 echo  $dir  $fq1 $fq2
 nohup $bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 & 
 cd ..
done
