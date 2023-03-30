cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc28/gatk_result/
for id in {PT1,PT2,PT4}
do 
grep -v 'Chr' ${id}.hg38_multianno.txt |cut -f 1-20|awk '{print $0"\t""'${id}'"}'  >${id}.maf ;
done



