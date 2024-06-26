cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/sh_cell_line/wgbs_20240611/methy_mtx



#bedtools intersect -a hg38_Huh7_H3K9ac.bedGraph -b dmr_sig.bedGraph -wo > dmr_sig_H3K9ac.bedGraph

#bedtools intersect -a hg38_Huh7_H3K4me1.bedGraph -b dmr_sig.bedGraph -wo > dmr_sig_H3K4me1.bedGraph

#bedtools intersect -a hg38_Huh7_H3K27ac.bedGraph -b dmr_sig.bedGraph -wo > dmr_sig_H3K27ac.bedGraph

#bedtools intersect -a hg38_Huh7_enhancer.bed -b dmr_sig.bedGraph -wo > dmr_sig_enhancer.bedGraph

bedtools intersect -a hg38_Huh7_H3K36ac.bedGraph -b dmr_sig.bedGraph -wo > dmr_sig_H3K36ac.bedGraph

bedtools intersect -a hg38_Huh7_H4K5ac.bedGraph -b dmr_sig.bedGraph -wo > dmr_sig_H4K5ac.bedGraph
