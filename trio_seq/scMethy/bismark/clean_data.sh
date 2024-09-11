cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/trio_seq/scMethy/raw_data/scMethy/hcc3/clean_data
for file in *.clean.fq.gz; do  
    # 使用正则表达式提取样本名（不包括 _1 或 _2 的部分）  
    if [[ $file =~ (.*_D[0-9]+)_[0-9]\.clean\.fq\.gz ]]; then  
        sample_name=${BASH_REMATCH[1]}  
        
        # 创建样本名对应的文件夹（如果不存在）  
        mkdir -p "$sample_name"  
        
        # 移动文件到对应的文件夹  
        mv "$file" "$sample_name/"  
    fi  
done 