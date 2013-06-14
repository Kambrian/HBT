#!/bin/bash
for run in 8213 8213.2 8213.4 8213.6 8213.8 8213.10 8213.12 8213.16
do
	make massfun_plot RUN_NUM=$run
	for i in 59 53 47 41 35
	do
		./massfun_plot $i 0
	done
done

