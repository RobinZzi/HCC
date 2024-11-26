cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/omix
for file in *.gz; do  
  #archive="${file%/}.tar.gz"  
  #tar -czf "$archive" "$dir" && \
  md5sum "$file" > "$file.md5" 
  #rm -rf "$dir"  
done