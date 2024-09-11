cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/trio_seq/bismark/merged_bam


for bam_file in *.bam; do  
    sorted_bam_file="${bam_file%.bam}_sorted.bam"  
    echo "Sorting $bam_file into $sorted_bam_file"  
    samtools sort -o "$sorted_bam_file" "$bam_file"  
    echo "Indexing $sorted_bam_file"  
    samtools index "$sorted_bam_file"  
done  