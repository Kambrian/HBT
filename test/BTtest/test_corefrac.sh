#!/bin/bash
snap=$1
fofid=$2
for CoreFrac in 0.01 0.04 0.111 0.25 0.51
do
./follow_fof.core $snap $fofid $CoreFrac &
done

for MassRelax in 10 5 3 2 1.4
do
	./follow_fof.adpt $snap $fofid $MassRelax 1 &
done
