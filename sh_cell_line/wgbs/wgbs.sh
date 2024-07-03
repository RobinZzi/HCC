folder=wgbs_20240701
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder
mkdir qc
mkdir trim_result
mkdir bismark_bam
mkdir remove_dup
mkdir meth_ex
mkdir cpg_report

cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/fastq
qcdir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/qc
path=$1
echo $path
files=$(ls $path)
for filename in $files
do
 cd $filename
 fastqc -t 12 -o $qcdir *.fastq.gz &
 cd ..
done

wait
echo "Current date and time: $(date)"
echo "All fastqc processes are completed."

cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/fastq
bin_trim_galore=trim_galore
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/trim_result/
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


wait
echo "Current date and time: $(date)"
echo "All trim processes are completed."



#bismark_bam
echo "PATH=/storage/zhangyanxiaoLab/zhangliwen/src/Bismark:$PATH" >> ${HOME}/.bashrc
source ~/.bashrc
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/trim_result
ref=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/bismark_bam
i=1
ls *fq.gz |while read id
    do
    headname=${id%_R*}
    if [ $i -eq 2 ];then
    echo ${headname}
    echo $i
    i=1
    bismark --genome $ref -1 ${headname}_R1_val_1.fq.gz -2 ${headname}_R2_val_2.fq.gz --o $dir &
    else
    i=$[$i+1]
    fi
    done

wait
echo "Current date and time: $(date)"
echo "All bismark_bam processes are completed."


#remove_dup
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/bismark_bam
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/remove_dup

cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
deduplicate_bismark --bam ${headname}.bam --output_dir $dir2 &
echo $headname
done


wait
echo "Current date and time: $(date)"
echo "All remove_dup processes are completed."




#meth_ex

dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/remove_dup
outP=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/meth_ex
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38


cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
bismark_methylation_extractor -p --no_overlap --parallel 30 --bedGraph --buffer_size 90G --cytosine_report --output ${outP} --genome_folder ${index} ${headname}.bam &
done



wait
echo "Current date and time: $(date)"
echo "All meth_ex processes are completed."

