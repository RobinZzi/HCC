cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/trim_result
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hisat2_genome/GRCh38_index
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/rna-seq/hisat2_result


i=1
ls *fq.gz |while read id
    do
    headname=${id%_R*}
    if [ $i -eq 2 ];then
    echo ${headname}
    echo $i
    i=1
    hisat2 -p 8 --dta -x $index -1 ${headname}_R1.fq.gz -2 ${headname}_R2.fq.gz -S ${dir}/${headname}.SAM &
    else
    i=$[$i+1]
    fi
    done
    
    