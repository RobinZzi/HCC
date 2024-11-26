cd /mnt/usb/hcc/HCC9/10x
for tar in *tar.gz; do
 md5sum "$tar" > "$tar.md5"
done

