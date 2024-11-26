cd /mnt/usb/hcc/HCC9/10x
#!/bin/bash  

# 创建一个空数组用于存储已处理的前缀  
processed_prefixes=()  

# 遍历所有 .fastq.gz 文件  
for file in *.fastq.gz; do  
  # 使用正则表达式提取前缀部分 HCC02-10X-PT5，以第一个下划线之前的部分为标准  
  prefix=$(echo "$file" | sed 's/^\(.*\)-[0-9]*_S[0-9]*_L[0-9]*_R[12]_.*\.fastq\.gz/\1/')  

  # 将提取的前缀由 '-' 替换为 '_'  
  new_prefix=$(echo "$prefix" | sed 's/-/_/g')  

  # 构造压缩文件的名字  
  archive_name="${new_prefix}_SCRNA_10x.tar.gz"  

  # 检查这个前缀是否已经被处理过  
  if [[ ! " ${processed_prefixes[@]} " =~ " ${prefix} " ]]; then  
    # 如果没有被处理过，将该前缀添加到已处理列表  
    processed_prefixes+=("$prefix")  

    # 找到所有该前缀相关的文件并打包成 tar.gz
    (tar -czf "$archive_name" "${prefix}-"*".fastq.gz"  
    echo "Created archive: $archive_name")& 
    echo "Creating archive: $archive_name" 
  fi  
done