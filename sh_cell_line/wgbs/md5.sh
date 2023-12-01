cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/sh_cell_line/wgbs
path=$1
files=$(ls$path)
for filename in $files
do
 cd $filename
 md5sum -c *.md5
 cd ..
done
