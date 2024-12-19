cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/hcc_multiomics_gsa_upload/data/WES/HCC2_WES

for file in *.gz;do
 (id="${file%.gz}"  
 gzip -d "$file"
 dos2unix "$id"
 gzip "$id")&
done