cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/hcc_multiomics_gsa_upload/scRNA/strt_Seq/
for dir in */; do  
    # Check if the item is a directory  
    if [ -d "$dir" ]; then  
        # Remove the trailing slash and rename the directory  
        mv "$dir" "${dir%/}_SCRNA_strt_Seq"  
    fi  
done  