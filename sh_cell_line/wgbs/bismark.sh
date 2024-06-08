cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/trim_result
ref=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38
dir=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_202405/bismark_bam
i=1
ls *fq.gz |while read id
    do
    headname=${id%_R*}
    if [ $i -eq 2 ];then
    echo ${headname}
    echo $i
    i=1
    bismark --genome $ref -1 ${headname}_R1.fq.gz -2 ${headname}_R2.fq.gz --o $dir &
    else
    i=$[$i+1]
    fi
    done



