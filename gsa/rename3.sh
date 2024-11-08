#!/bin/bash  
cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/hcc_multiomics_gsa_upload/scPBAT/HCC6
# Iterate over all files that contain 'PTnt' in their name  
for file in *PTnt*; do  
    # Generate the new filename by replacing 'PTnt' with 'NT'  
    new_filename="${file//PTnt/NT}"  
    
    # Rename the file  
    if [ "$file" != "$new_filename" ]; then  
        echo "Renaming '$file' to '$new_filename'"  
        mv "$file" "$new_filename"  
    fi  
done