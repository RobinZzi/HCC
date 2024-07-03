folder=wgbs_20240701
sample_num=4
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder
mkdir qc
mkdir trim_result
mkdir bismark_bam
mkdir remove_dup
mkdir meth_ex
mkdir cpg_report



cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/fastq
bin_trim_galore=trim_galore
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/trim_result/
tmp_dir=$(mktemp -d)
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
 ($bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
  touch "$tmp_dir/$headname.done") &
 done
 cd ..
done

echo "#####################################"
echo "Current date and time: $(date)"
echo "All trim processes are in processing."
echo "#####################################"
       
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

echo "#####################################"
echo "Current date and time: $(date)"
echo "All trim processes are completed."
echo "#####################################"


#bismark_bam
echo "PATH=/storage/zhangyanxiaoLab/zhangliwen/src/Bismark:$PATH" >> ${HOME}/.bashrc
source ~/.bashrc
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/trim_result
ref=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/bismark_bam
tmp_dir=$(mktemp -d)

i=1
ls *fq.gz |while read id
    do
    headname=${id%_R*}
    if [ $i -eq 2 ];then
    echo ${headname}
    echo $i
    i=1
    (bismark --genome $ref -1 ${headname}_R1_val_1.fq.gz -2 ${headname}_R2_val_2.fq.gz --o $dir
    touch "$tmp_dir/$headname.done") &
    else
    i=$[$i+1]
    fi
    done

echo "#####################################"
echo "Current date and time: $(date)"
echo "All bismark_bam are in processing."
echo "#####################################"
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

echo "#####################################"
echo "Current date and time: $(date)"
echo "All bismark_bam processes are completed."
echo "#####################################"

#remove_dup
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/bismark_bam
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/remove_dup

cd $dir
tmp_dir=$(mktemp -d)
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
(deduplicate_bismark --bam ${headname}.bam --output_dir $dir2
 touch "$tmp_dir/$headname.done") &
echo $headname
done


echo "#####################################"
echo "Current date and time: $(date)"
echo "All remove_dup are in processing."
echo "#####################################"
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

wait
echo "#####################################"
echo "Current date and time: $(date)"
echo "All remove_dup processes are completed."
echo "#####################################"




#meth_ex

dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/remove_dup
outP=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/meth_ex
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38


cd $dir

ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
(bismark_methylation_extractor -p --no_overlap --parallel 30 --bedGraph --buffer_size 90G --cytosine_report --output ${outP} --genome_folder ${index} ${headname}.bam 
 touch "$tmp_dir/$headname.done") &
done


echo "#####################################"
echo "Current date and time: $(date)"
echo "All meth_ex are in processing."
echo "#####################################"

while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir
echo "#####################################"
echo "Current date and time: $(date)"
echo "All meth_ex processes are completed."
echo "#####################################"





cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/meth_ex
mkdir CpG_report
mkdir bismark_cov

mv *.CpG_report.txt CpG_report
mv *.bismark.cov.gz bismark_cov



