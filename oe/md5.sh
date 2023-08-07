cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/oe/ex3
path=$1
files=$(ls$path)
for filename in $files
do
 cd $filename
 md5sum -c *.md5
 cd ..
done