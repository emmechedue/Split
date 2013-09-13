#!/bin/bash

mkdir configs

for i in $(seq 1 1 121)
do
	J=$((i-1))
	esse=`echo $J \\* 0.05 | bc`
	sed "s|s= 0.05|s= $esse|" ./dummy.conf >./configs/$i.conf
done


