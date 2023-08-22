cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report/10kbin
bedtools intersect -a g6_10k_global_regions.bed -b H3K9me3.bed -wo >g6_global_H3K9me3.bedGraph
bedtools intersect -a g6_10k_global_regions.bed -b H3K27me3.bed -wo >g6_global_H3K27me3.bedGraph
bedtools intersect -a g6_10k_global_regions.bed -b H3K36me3.bed -wo >g6_global_H3K36me3.bedGraph
bedtools intersect -a g6_10k_global_regions.bed -b H3K27ac.bed -wo >g6_global_H3K27ac.bedGraph

bedtools intersect -a g6_10k_dmr_regions.bed -b H3K9me3.bed -wo >g6_dmr_H3K9me3.bedGraph
bedtools intersect -a g6_10k_dmr_regions.bed -b H3K27me3.bed -wo >g6_dmr_H3K27me3.bedGraph
bedtools intersect -a g6_10k_dmr_regions.bed -b H3K36me3.bed -wo >g6_dmr_H3K36me3.bedGraph
bedtools intersect -a g6_10k_dmr_regions.bed -b H3K27ac.bed -wo >g6_dmr_H3K27ac.bedGraph 

bedtools intersect -a g6_10k_admr_regions.bed  -b H3K9me3.bed -wo >g6_admr_H3K9me3.bedGraph
bedtools intersect -a g6_10k_admr_regions.bed -b H3K27me3.bed -wo >g6_admr_H3K27me3.bedGraph
bedtools intersect -a g6_10k_admr_regions.bed -b H3K36me3.bed -wo >g6_admr_H3K36me3.bedGraph
bedtools intersect -a g6_10k_admr_regions.bed -b H3K27ac.bed -wo >g6_admr_H3K27ac.bedGraph

bedtools intersect -a g6_10k_admr_regions.bed -b PMD_coordinates_hg38.bed  -wo > g6_admr_PMD.bedGraph
bedtools intersect -a g6_10k_dmr_regions.bed -b PMD_coordinates_hg38.bed  -wo > g6_dmr_PMD.bedGraph
bedtools intersect -a g6_10k_global_regions.bed -b PMD_coordinates_hg38.bed  -wo > g6_global_PMD.bedGraph