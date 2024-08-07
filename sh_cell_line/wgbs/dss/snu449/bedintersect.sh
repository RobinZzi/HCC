genome_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/liver_cancer_cellline/chip_seq/HepG2
data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs/dss/SNU449/dmr_anno
data=GA45_dmrs.bedGraph
cd $data_path






bedtools intersect -a $genome_path/PMD_coordinates_hg38.bed -b $data -wo > dmr_pmd.bedGraph

bedtools intersect -a $genome_path/hg38_H3K4me1.bed.gz -b $data -wo > dmr_H3K4me1.bedGraph

bedtools intersect -a $genome_path/hg38_H3K9me3.bed.gz -b $data -wo > dmr_H3K9me3.bedGraph

bedtools intersect -a $genome_path/hg38_H3K36me3.bed.gz -b $data -wo > dmr_H3K36me3.bedGraph

bedtools intersect -a $genome_path/hg38_H3K27ac.bed.gz -b $data -wo > dmr_H3K27ac.bedGraph

bedtools intersect -a $genome_path/hg38_H3K9ac.bed.gz -b $data -wo > dmr_H3K9ac.bedGraph

bedtools intersect -a $genome_path/hg38_H3K4me3.bed.gz -b $data -wo > dmr_H3K4me3.bedGraph

bedtools intersect -a $genome_path/hg38_H3K27me3.bed.gz -b $data -wo > dmr_H3K27me3.bedGraph

bedtools intersect -a $genome_path/hg38_promoter.bed -b $data -wo > dmr_promoter.bedGraph





