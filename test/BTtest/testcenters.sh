#!/bin/bash
# 8213:
#	S40G2017 sub31462
#	S37G2115 sub28224
# 6702DM: S51G86 sub7990
#  	  S43G11 sub5562

snap=$1
fofid=$2
echo snap=$1 fofid=$2

echo mean
./follow_fof.mean $1 $2 
echo minpotH:
./follow_fof.minpotH $1 $2 
echo potW:
./follow_fof.potW $1 $2 

for CoreFrac in 0.01 0.04 0.111 0.25 0.51
do
echo core$CoreFrac$
./follow_fof.core $snap $fofid $CoreFrac 
done

