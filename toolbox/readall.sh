#!/bin/bash

for ((i=0;i<100;i++))
do
	if ./loadview $i 1>>6402load.log 2>>6402err.log
	then
		echo $i
	else
		echo wrong $i
	fi	
done	
