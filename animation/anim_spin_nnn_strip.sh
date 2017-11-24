#!/bin/bash


#num = 5
declare -r num=10 
declare -r kpt=1000
declare -r steps=100

#echo "$num"

#VAR=$(echo "scale=2; $IMG_WIDTH/$IMG2_WIDTH" | bc)

for i in $(seq 0 $steps);
do
  #echo "${i}"
  VAR=$(bc <<<"scale=4; ${i}/(20*$steps)")
  VAR2=$(($steps - ${i}))
  #VAR=$(bc -l <<<"${i}/$steps")
  #echo "$VAR"
  ./tight_binding.x -n $num -k $kpt -edge 1.0 -t2 $VAR -s 1.0 > u${i}.dat
  ./tight_binding.x -n $num -k $kpt -edge 1.0 -t2 $VAR -s -1.0 > d${i}.dat
  python band_plot_spin.py $num $kpt u${i}.dat d${i}.dat out${i}.jpg
done

for i in $(seq 1 $steps);
do
  #echo "${i}"
  VAR=$(bc <<<"scale=4; 1-(${i}/$steps)")
  #echo "$VAR"
  #VAR2=$(($steps - ${i}))
  VAR2=$((${i} + $steps))
  #VAR=$(bc -l <<<"${i}/$steps")
  #echo "$VAR2"
  ./tight_binding.x -n $num -k $kpt -edge $VAR -t2 0.05 -s 1.0> u$VAR2.dat
  ./tight_binding.x -n $num -k $kpt -edge $VAR -t2 0.05 -s -1.0> d$VAR2.dat
  python band_plot_spin.py $num $kpt u$VAR2.dat d$VAR2.dat out$VAR2.jpg
done


ffmpeg -i out%d.jpg -vcodec libx264 -s 640x480 -pix_fmt yuv420p turnon_nnn_0.05_turn_on_strip_withspin.mp4
rm *.dat
rm *.jpg