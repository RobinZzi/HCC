cd /storage/zhangyanxiaoLab/zhangliwen/projects/hcc/data/oe
path=$1
files=$(ls$path)
for filename in $files
do
 md5sum -c *.md5
 cd ..
done