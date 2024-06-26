folder=wgbs_20240614
cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder
mkdir qc
mkdir trim_result
mkdir bismark_bam
mkdir remove_dup
mkdir meth_ex
mkdir cpg_report



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
    bismark --genome $ref -1 ${headname}_R1_val_1.fq.gz -2 ${headname}_R2_val_2.fq.gz --o $dir 
    else
    i=$[$i+1]
    fi
    done

wait


#remove_dup
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/bismark_bam
dir2=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/remove_dup

cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
deduplicate_bismark --bam ${headname}.bam --output_dir $dir2
echo $headname
done


wait




#meth_ex

dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/remove_dup
outP=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/$folder/meth_ex
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38


cd $dir
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
bismark_methylation_extractor -p --no_overlap --parallel 30 --bedGraph --buffer_size 90G --cytosine_report --output ${outP} --genome_folder ${index} ${headname}.bam
done

x



