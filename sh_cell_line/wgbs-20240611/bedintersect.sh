genome_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/liver_cancer_cellline/chip_seq/Huh7
data_path=/storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240611/methy_mtx/dmr_anno/
data=SNG6_dmrs.bedGraph
cd $data_path




bedtools intersect -a $genome_path/PMD_coordinates_hg38.bed -b $data -wo > dmr_pmd.bedGraph

bedtools intersect -a $genome_path/hg38_Huh7_H4K5ac.bedGraph -b $data -wo > dmr_H4K5ac.bedGraph

bedtools intersect -a $genome_path/hg38_Huh7_H3K9ac.bedGraph -b $data -wo > dmr_H3K9ac.bedGraph

bedtools intersect -a $genome_path/hg38_Huh7_H3K27ac.bedGraph -b $data -wo > dmr_H3K27ac.bedGraph

bedtools intersect -a $genome_path/hg38_Huh7_H3K36ac.bedGraph -b $data -wo > dmr_H3K36ac.bedGraph

bedtools intersect -a $genome_path/hg38_Huh7_H3K4me1.bedGraph -b $data -wo > dmr_H3K4me1.bedGraph

bedtools intersect -a $genome_path/hg38_Huh7_enhancer.bed -b $data -wo > dmr_enhancer.bedGraph

bedtools intersect -a $genome_path/hg38_promoter.bed -b $data -wo > dmr_promoter.bedGraph




