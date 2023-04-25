cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/TCGA_related/methy_merge
bedtools intersect -a H3K9me3.bed -b late_demeth_regions.bedGraph -wo > late_H3K9me3.bedGraph
bedtools intersect -a H3K27me3.bed -b late_demeth_regions.bedGraph -wo > late_H3K27me3.bedGraph
bedtools intersect -a H3K4me3.bed -b late_demeth_regions.bedGraph -wo > late_H3K4me3.bedGraph
bedtools intersect -a H3K36me3.bed -b late_demeth_regions.bedGraph -wo > late_H3K36me3.bedGraph
bedtools intersect -a H3K27ac.bed -b late_demeth_regions.bedGraph -wo > late_H3K27ac.bedGraph



bedtools intersect -a H3K9me3.bed -b early_demeth_regions.bedGraph -wo > early_H3K9me3.bedGraph
bedtools intersect -a H3K27me3.bed -b early_demeth_regions.bedGraph -wo > early_H3K27me3.bedGraph
bedtools intersect -a H3K4me3.bed -b early_demeth_regions.bedGraph -wo > early_H3K4me3.bedGraph
bedtools intersect -a H3K36me3.bed -b early_demeth_regions.bedGraph -wo > early_H3K36me3.bedGraph
bedtools intersect -a H3K27ac.bed -b early_demeth_regions.bedGraph -wo > early_H3K27ac.bedGraph