cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240701/bismark_cov



bedtools intersect -a H3K9me3.bed.gz -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K9me3.bedGraph

bedtools intersect -a H3K4me3.bed.gz -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K4me3.bedGraph

bedtools intersect -a H3K27ac.bed.gz -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K27ac.bedGraph

bedtools intersect -a H3K36me3.bed.gz -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K36me3.bedGraph

bedtools intersect -a hg38_promoter.bed -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_promoter.bedGraph

bedtools intersect -a PMD_coordinates_hg38.bed -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_pmd.bedGraph

bedtools intersect -a hg38_293T_H3K27ac.bedGraph -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K27ac_2.bedGraph

bedtools intersect -a hg38_293T_H3K4me3.bedGraph -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K4me3_2.bedGraph

bedtools intersect -a 293T_H3K4me3.bedgraph -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K4me3_3.bedGraph

bedtools intersect -a 293T_H3K9me3.bedgraph -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K9me3_2.bedGraph



bedtools intersect -a H3K9me3_3.bed.gz -b GA45_dmrs_sig.bedGraph -wo > dmr_sig_H3K9me3_3.bedGraph