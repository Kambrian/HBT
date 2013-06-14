#!/bin/bash
snap=$1
fofid=$2
for skip in 2 5 10
do
	./follow_fof.adpt $snap $fofid 2 $skip &
done
