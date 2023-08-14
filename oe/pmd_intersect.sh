cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/analysis/oe/merge2/methy_report

bedtools intersect -a ga_admr_regions.bedGraph -b PMD_coordinates_hg38.bed  -wo > ga_admr_PMD.bedGraph
bedtools intersect -a ga_dmr_regions.bedGraph -b PMD_coordinates_hg38.bed  -wo > ga_dmr_PMD.bedGraph
bedtools intersect -a ga_global_regions.bedGraph -b PMD_coordinates_hg38.bed  -wo > ga_global_PMD.bedGraph