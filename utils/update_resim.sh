#!/bin/bash
for i in 6[4-7]0[0-9]
#for i in 6501 6506 6600 6601 6602 6700 6701 6702
do	
echo $i start
cd ~/data/$i/source/v8.3
chmod u+x autorun
./autorun $i 1>autorun.$i.log 2>&1
date
echo $i done
done

