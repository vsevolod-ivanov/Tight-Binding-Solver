#!/bin/bash


#num = 5
declare -r num=10 
declare -r kpt=100
declare -r steps=10

#echo "$num"

#VAR=$(echo "scale=2; $IMG_WIDTH/$IMG2_WIDTH" | bc)

for i in $(seq 0 $steps);
do
  #echo "${i}"
  VAR=$(bc <<<"scale=2; ${i}/$steps")
  VAR2=$(($steps - ${i}))
  #VAR=$(bc -l <<<"${i}/$steps")
  #echo "$VAR"
  ./tight_binding.x -n $num -k $kpt -edge $VAR > ${i}.dat
  python band_plot.py $num $kpt ${i}.dat out$VAR2.jpg
done
