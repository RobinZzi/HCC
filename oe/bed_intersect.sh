cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge
bedtools intersect -a ga_dmr_regions.bedGraph -b H3K9me3.bed -wo >dmr_H3K9me3.bedGraph
bedtools intersect -a ga_dmr_regions.bedGraph -b H3K27me3.bed -wo >dmr_H3K27me3.bedGraph
bedtools intersect -a ga_dmr_regions.bedGraph -b H3K36me3.bed -wo >dmr_H3K36me3.bedGraph
bedtools intersect -a ga_dmr_regions.bedGraph -b H3K27ac.bed -wo >dmr_H3K27ac.bedGraph 