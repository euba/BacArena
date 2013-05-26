mkdir temp
cp *.jpeg temp/.
mogrify -resize 800x800  temp/*.JPG
convert temp/*.jpeg -delay 10 -morph 10 temp/%05d.jpeg
ffmpeg -r 25 -qscale 2  -i temp/%05d.jpeg output.mp4
rm -R temp


