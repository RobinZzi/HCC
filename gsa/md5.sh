(cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/GSA_backup/hcc_multiomics_gsa_upload/CellLine/OverExpress/WGBS/HepG2_OE_SNHG6
for tar in *.gz; do
 md5sum "$tar" > "$tar.md5"
done)&

(cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/GSA_backup/hcc_multiomics_gsa_upload/CellLine/OverExpress/WGBS/HepG2_OE_GADD45A
for tar in *.gz; do
 md5sum "$tar" > "$tar.md5"
done)&

(cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/GSA_backup/hcc_multiomics_gsa_upload/CellLine/KnockDown/WGBS/Huh7_KD_SNHG6_WGBS
for tar in *.gz; do
 md5sum "$tar" > "$tar.md5"
done)&

(cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/GSA_backup/hcc_multiomics_gsa_upload/CellLine/KnockDown/WGBS/293T_KD_GADD45A_WGBS
for tar in *.gz; do
 md5sum "$tar" > "$tar.md5"
done)&

(cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/GSA_backup/hcc_multiomics_gsa_upload/CellLine/KnockDown/RNA_Seq/Huh7_KD_SNHG6_RNA_Seq
for tar in *.gz; do
 md5sum "$tar" > "$tar.md5"
done)&

(cd /mnt/transposon1/zhangyanxiaoLab/zhangliwen/GSA_backup/hcc_multiomics_gsa_upload/CellLine/KnockDown/RNA_Seq/293T_KD_GADD45A_RNA_Seq
for tar in *.gz; do
 md5sum "$tar" > "$tar.md5"
done)&