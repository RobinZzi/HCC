folder=hepg2_trioseq_2
data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/cell_re_rrbs/$folder
sample_num=6
cd $data_path
mkdir qc
mkdir trim_result
mkdir bismark_bam
mkdir remove_dup
mkdir meth_ex
mkdir cpg_report

cd $data_path/fastq
qcdir=$data_path/qc
tmp_dir=$(mktemp -d)

path=$1
echo $path
files=$(ls $path)
for filename in $files
do
 cd $filename
 (fastqc -t 12 -o $qcdir *.fastq.gz 
 touch "$tmp_dir/$filename.done") &
 cd ..
done

echo "#####################################"
echo "Current date and time: $(date)"
echo "All qc processes are in processing."
echo "#####################################"
       
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

echo "#####################################"
echo "Current date and time: $(date)"
echo "All qc processes are completed."
echo "#####################################"

cd $data_path/fastq
bin_trim_galore=trim_galore
dir=$data_path/trim_result/
tmp_dir=$(mktemp -d)
path=$1
files=$(ls$path)

for filename in $files
do
 cd $filename
 ls *fastq.gz|while read id
 do
 fq1=${filename}_R1.fastq.gz
 fq2=${filename}_R2.fastq.gz
 echo  $fq1 $fq2
 ($bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
  touch "$tmp_dir/$filename.done") &
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
cd $data_path/trim_result
ref=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
dir=$data_path/bismark_bam
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
dir=$data_path/bismark_bam
dir2=$data_path/remove_dup

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

dir=$data_path/remove_dup
outP=$data_path/meth_ex
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38


cd $dir
tmp_dir=$(mktemp -d)
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





cd $data_path/meth_ex
mkdir CpG_report
mkdir bismark_cov

mv *.CpG_report.txt CpG_report
mv *.bismark.cov.gz bismark_cov



