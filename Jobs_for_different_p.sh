#!/bin/bash
PMAX=10
PMIN=0
PSTEP=1

for i in $(seq $PMIN $PSTEP $PMAX)
do
	mkdir -p p+$i
	sed "s|p= 10.0|p= $i.0|" ./dummy.conf >./p+$i/dummy.conf
	cp create_configs.sh ./p+$i/create_configs.sh
	cp multiple.sge ./p+$i/multiple.sge
	cp a.out ./p+$i/a.out
	cd ./p+$i
	sh create_configs.sh
	cd ..
done

