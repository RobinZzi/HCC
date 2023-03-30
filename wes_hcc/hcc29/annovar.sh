cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/wes_hcc/hcc29/gatk_result/
for id in {PT1,PT3,PT4}
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