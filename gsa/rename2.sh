#!/bin/bash  
cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/hcc_multiomics_gsa_upload/scPBAT/HCC6
# Iterate over all .clean.fq.gz files  
for file in *.clean.fq.gz; do  
    # Extract parts of the filename using bash pattern matching  
    # Assume filenames are like pt1_D194_1.clean.fq.gz  
    pt_number=${file%%_*}  # pt1  
    rest_of_file=${file#*_}  # D194_1.clean.fq.gz  
    d_number=${rest_of_file%%_*}  # D194  
    read_number_with_extension=${rest_of_file##*_}  # 1.clean.fq.gz  
    read_number=${read_number_with_extension%%.*}  # 1  

    # Construct the new filename using extracted parts  
    # Modify "HCC?" as needed  
    new_filename="HCC6_PT${pt_number#pt}_D${d_number#D}_SCPBAT_R${read_number}.fq.gz"  

    # Rename the file  
    echo "Renaming '$file' to '$new_filename'"  
    mv "$file" "$new_filename"  
done