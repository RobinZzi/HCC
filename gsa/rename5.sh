#!/bin/bash  

# 根目录  
base_dir="/mnt/usb/hcc/HCC2/strt"  

# 找到所有符合命名模式的文件并处理  
find "$base_dir" -type f -name "*.fastq.gz" | while read -r file; do  
# 提取目录名以便构建新文件名  
# 比如：从 /mnt/usb/hcc/HCC3/strt 中提取 HCC3 和 strt  
dir=$(dirname "$file")  
hcc_name=$(basename "$(dirname "$dir")")  
drop_seq=$(basename "$dir")  

# 提取文件名中的信息  
base_name=$(basename "$file" .fastq.gz)  
pt_number=$(echo "$base_name" | cut -d'_' -f1 | sed 's/pt/PT/')  

# 如果样本编号以 'mt' 开头，则替换为 'NT'  
pt_number=$(echo "$pt_number" | sed 's/^MT/NT/')  

# 根据 pt1_ 后的数字决定 R1 或 R2  
number_suffix=$(echo "$base_name" | cut -d'_' -f2)  
if [[ $number_suffix -eq 1 ]]; then  
r_suffix="R1"  
elif [[ $number_suffix -eq 2 ]]; then  
r_suffix="R2"  
else  
  echo "Unknown suffix for file: $file"  
continue  
fi  

# 创建新的文件名  
new_name="${hcc_name}_${pt_number}_SCRNA_${drop_seq}_${r_suffix}.fq.gz"  

# 重命名文件  
mv "$file" "$dir/$new_name"  
echo "Renamed: $file to $dir/$new_name"  
done