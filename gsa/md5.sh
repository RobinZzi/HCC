cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/hcc_multiomics_gsa_upload/data/WES/HCC9_WES
for file in *.gz; do
 md5sum "$file" > "$file.md5" &
done

