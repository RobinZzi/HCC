cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/fastq
bin_trim_galore=trim_galore
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/trim_result/
path=$1
files=$(ls$path)
for filename in $files
do
 cd $filename
 ls *fastq.gz|while read id
 do
 headname=${id%_R*}
 fq1=${headname}_R1.fastq.gz
 fq2=${headname}_R2.fastq.gz
 echo  $fq1 $fq2
 nohup $bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 &
 done
 cd ..
done
