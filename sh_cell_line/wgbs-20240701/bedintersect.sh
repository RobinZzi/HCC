cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_cov



bedtools intersect -a H3K9me3.bed.gz -b GA45_dmrs.bedGraph -wo > dmr_H3K9me3.bedGraph

bedtools intersect -a H3K4me3.bed.gz -b GA45_dmrs.bedGraph -wo > dmr_H3K4me3.bedGraph

bedtools intersect -a H3K27ac.bed.gz -b GA45_dmrs.bedGraph -wo > dmr_H3K27ac.bedGraph

bedtools intersect -a H3K36me3.bed.gz -b GA45_dmrs.bedGraph -wo > dmr_H3K36me3.bedGraph

bedtools intersect -a hg38_promoter.bed -b GA45_dmrs.bedGraph -wo > dmr_promoter.bedGraph

bedtools intersect -a PMD_coordinates_hg38.bed -b GA45_dmrs.bedGraph -wo > dmr_pmd.bedGraph
