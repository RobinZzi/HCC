bin_trim_galore=trim_galore
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/over_expression/trim_result/
i=1
ls *fastq.gz|while read id
do
if [ $i -eq 2 ];then
    i=1
else
headname=${id%.R*}
fq1=${headname}.R1.fastq.gz
fq2=${headname}.R2.fastq.gz
echo  $dir  $fq1 $fq2
nohup $bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 & 
i=$[$i+1]
fi
done 

