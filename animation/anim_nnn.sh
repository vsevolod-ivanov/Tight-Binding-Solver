#!/bin/bash


#num = 5
declare -r num=10 
declare -r kpt=1000
declare -r steps=5

#echo "$num"

#VAR=$(echo "scale=2; $IMG_WIDTH/$IMG2_WIDTH" | bc)

for i in $(seq 0 $steps);
do
  #echo "${i}"
  VAR=$(bc <<<"scale=4; ${i}/(20*$steps)")
  VAR2=$(($steps - ${i}))
  #VAR=$(bc -l <<<"${i}/$steps")
  #echo "$VAR"
  ./tight_binding.x -n $num -k $kpt -edge 1.0 -t2 $VAR > ${i}.dat
  python band_plot.py $num $kpt ${i}.dat out${i}.jpg
done

ffmpeg -i out%d.jpg -vcodec libx264 -s 640x480 -pix_fmt yuv420p turnon_nnn_0.05_nostrip.mp4
#rm *.dat
#rm *.jpg