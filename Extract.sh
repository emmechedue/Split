#!/bin/bash

NUMJOB=36278   # Here I have to give the script the number of the job!

for i in $(seq 1 1 151)
do
	j=`echo $i-1 | bc`
	esse=`echo $j \\* 0.01 | bc`		#Computing the value of s
	dummy0=$(echo "$esse == 0" | bc)
	dummy=$(echo "$esse < 1" | bc)		#Here I create the folder named s+... or s+0....
	if [ $dummy0 -eq 1 ]; then
		DIR=s+0.00
	elif [ $dummy -eq 1 ]; then
		DIR=s+0$esse
	else
		DIR=s+$esse
	fi
	mkdir -p $DIR
	tar -xf job$NUMJOB.$i.tgz -C ./$DIR		#Unzipping the file
	mv ./$DIR/data/Stefano.Duca/job$NUMJOB.$i/ensambleN.txt ./$DIR/data/Stefano.Duca/job$NUMJOB.$i/ensamblex.txt ./$DIR/data/Stefano.Duca/job$NUMJOB.$i/parameters.txt ./$DIR/data/Stefano.Duca/job$NUMJOB.$i/time.txt ./$DIR		#Putting the files in the proper folder and deleting everything that I don't need
	rm -rf ./$DIR/data	
	rm ./job$NUMJOB.$i.tgz -f
done

python Multiple_analysis.py

#find -name *ensambleN.txt* -delete
#find -name *ensamblex.txt* -delete

