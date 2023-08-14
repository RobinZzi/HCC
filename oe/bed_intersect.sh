cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report
bedtools intersect -a ga_admr_regions.bedGraph -b H3K9me3.bed -wo >ga_admr_H3K9me3.bedGraph
bedtools intersect -a ga_admr_regions.bedGraph -b H3K27me3.bed -wo >ga_admr_H3K27me3.bedGraph
bedtools intersect -a ga_admr_regions.bedGraph -b H3K36me3.bed -wo >ga_admr_H3K36me3.bedGraph
bedtools intersect -a ga_admr_regions.bedGraph -b H3K27ac.bed -wo >ga_admr_H3K27ac.bedGraph

bedtools intersect -a g6_admr_regions.bedGraph -b H3K9me3.bed -wo >g6_admr_H3K9me3.bedGraph
bedtools intersect -a g6_admr_regions.bedGraph -b H3K27me3.bed -wo >g6_admr_H3K27me3.bedGraph
bedtools intersect -a g6_admr_regions.bedGraph -b H3K36me3.bed -wo >g6_admr_H3K36me3.bedGraph
bedtools intersect -a g6_admr_regions.bedGraph -b H3K27ac.bed -wo >g6_admr_H3K27ac.bedGraph 

bedtools intersect -a ga_dmr_regions.bedGraph -b H3K9me3.bed -wo >ga_dmr_H3K9me3.bedGraph
bedtools intersect -a ga_dmr_regions.bedGraph -b H3K27me3.bed -wo >ga_dmr_H3K27me3.bedGraph
bedtools intersect -a ga_dmr_regions.bedGraph -b H3K36me3.bed -wo >ga_dmr_H3K36me3.bedGraph
bedtools intersect -a ga_dmr_regions.bedGraph -b H3K27ac.bed -wo >ga_dmr_H3K27ac.bedGraph

bedtools intersect -a g6_dmr_regions.bedGraph -b H3K9me3.bed -wo >g6_dmr_H3K9me3.bedGraph
bedtools intersect -a g6_dmr_regions.bedGraph -b H3K27me3.bed -wo >g6_dmr_H3K27me3.bedGraph
bedtools intersect -a g6_dmr_regions.bedGraph -b H3K36me3.bed -wo >g6_dmr_H3K36me3.bedGraph
bedtools intersect -a g6_dmr_regions.bedGraph -b H3K27ac.bed -wo >g6_dmr_H3K27ac.bedGraph 

bedtools intersect -a ga_global_regions.bedGraph -b H3K9me3.bed -wo >ga_global_H3K9me3.bedGraph
bedtools intersect -a ga_global_regions.bedGraph -b H3K27me3.bed -wo >ga_global_H3K27me3.bedGraph
bedtools intersect -a ga_global_regions.bedGraph -b H3K36me3.bed -wo >ga_global_H3K36me3.bedGraph
bedtools intersect -a ga_global_regions.bedGraph -b H3K27ac.bed -wo >ga_global_H3K27ac.bedGraph

bedtools intersect -a g6_global_regions.bedGraph -b H3K9me3.bed -wo >g6_global_H3K9me3.bedGraph
bedtools intersect -a g6_global_regions.bedGraph -b H3K27me3.bed -wo >g6_global_H3K27me3.bedGraph
bedtools intersect -a g6_global_regions.bedGraph -b H3K36me3.bed -wo >g6_global_H3K36me3.bedGraph
bedtools intersect -a g6_global_regions.bedGraph -b H3K27ac.bed -wo >g6_global_H3K27ac.bedGraph 

