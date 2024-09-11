
sample_num=46
folder=pt2
#bismark_bam

cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/trim_result/${folder}
ref=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/bismark_bam/${folder}
tmp_dir=$(mktemp -d)

i=1
ls *fq.gz |while read id
    do
    headname=${id%_R*}
    if [ $i -eq 2 ];then
    echo ${headname}
    echo $i
    i=1
    (bismark --genome $ref -1 ${headname}_R1.fq.gz -2 ${headname}_R2.fq.gz --o $dir --pbat
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
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/bismark_bam/${folder}
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/remove_dup/${folder}

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






