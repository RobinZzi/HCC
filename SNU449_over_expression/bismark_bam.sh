ref=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/over_expression/bismark_bam
i=1
ls *fq.gz |while read id
    do
    headname=${id%.R*}
    if [ $i -eq 2 ];then
    echo ${headname}
    echo $i
    i=1
    bismark --genome $ref -1 ${headname}.R1_val_1.fq.gz -2 ${headname}.R2_val_2.fq.gz --o $dir &
    else
    i=$[$i+1]
    fi
    done



