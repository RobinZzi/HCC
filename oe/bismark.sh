cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/trim_result
ref=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/bismark_bam
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



