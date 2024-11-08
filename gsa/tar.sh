cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/GSA_backup/hcc_multiomics_gsa_upload/scRNA/strt_Seq
for dir in */; do  
  archive="${dir%/}.tar.gz"  
  tar -czf "$archive" "$dir" && \
  md5sum "$archive" > "$archive.md5" && \
  rm -rf "$dir"  
done