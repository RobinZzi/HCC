cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/wes_hcc/hcc3/gatk_result/bqsr_result

for id in {PT1,PT2,PT3,PT4}
do
/storage/zhangyanxiaoLab/zhangliwen/src/annovar/table_annovar.pl ${id}_filter.vcf /storage/zhangyanxiaoLab/zhangliwen/genome/humandb/ \
-buildver hg38 \
-out ${id} \
-remove \
-protocol refGene,knownGene,clinvar_20170905 \
-operation g,g,f \
-nastring . \
-vcfinput
done