#!/bin/bash
rm -f msplot.log

for((i=0;i<100;i+=1))
do
for((j=0;j<3;j+=1))
do	
echo "massfun_plot $i $j:" >>msplot.log 2>&1
./massfun_plot $i $j 1>>msplot.log 2>&1
done
done
