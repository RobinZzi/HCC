genome_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240705/bismark_cov/genome_bed
data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240705/bismark_cov/dmr_dmcp0.01delta0.2
data=GA45_dmrs.bedGraph
cd $data_path




bedtools intersect -a $genome_path/PMD_coordinates_hg38.bed -b $data -wo > dmr_pmd.bedGraph

bedtools intersect -a $genome_path/hg38_293T_H3K9me3.bedGraph -b $data -wo > dmr_H3K9me3.bedGraph

bedtools intersect -a $genome_path/hg38_293T_H3K4me3.bedGraph -b $data -wo > dmr_H3K4me3.bedGraph

bedtools intersect -a $genome_path/hg38_293T_H3K27ac.bedGraph -b $data -wo > dmr_H3K27ac.bedGraph

bedtools intersect -a $genome_path/hg38_293T_H3K4me3_2.bedGraph -b $data -wo > dmr_H3K4me3_2.bedGraph

bedtools intersect -a $genome_path/H3K27ac.bed.gz -b $data -wo > dmr_H3K27ac_2.bedGraph

bedtools intersect -a $genome_path/H3K9me3.bed.gz -b $data -wo > dmr_H3K9me3_2.bedGraph

bedtools intersect -a $genome_path/H3K36me3.bed.gz -b $data -wo > dmr_H3K36me3.bedGraph




