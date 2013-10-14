#!/bin/sh

i=0
while [ $i -le 7 ]; do

  ./plot_all5.sh $i
  ./plot_all.sh $i
  ./plot_all3.sh $i
  ./plot_all2.sh $i
  ./plot_all4.sh $i
  ./plot_all6.sh $i

  i=$((i+1))
done

exit
