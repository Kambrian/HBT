#!/bin/bash

export OMP_NUM_THREADS=50

for run in 6114  6115  6116  6117  6119  6120  6121  6122  6123  6127  6128  6129  6140  6141  6142  6143  6144
do
	for((i=0;i<100;i+=1))
	do
	echo "Snap=$i:" >>log/makeprof.$run.log 2>&1
	./haloprof.$run $i 1>>log/makeprof.$run.log 2>&1
	done
	chmod a-w ~/data/$run/subcat/profile/logbin/*
done
