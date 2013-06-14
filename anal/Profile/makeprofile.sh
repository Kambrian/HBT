#!/bin/bash
rm -f makeprof.log
rm -f makedynprof.sph.log

for((i=0;i<100;i+=1))
do
echo "haloprof $i:" >>makeprof.log 2>&1
../haloprof $i 1>>makeprof.log 2>&1
echo "gashaloprof $i:" >>makeprof.log 2>&1
../gashaloprof $i >>makeprof.log 2>&1
echo "dynamhaloprof.sph $i:" >> makedynprof.sph.log 2>&1
../dynamprof.sph $i >> makedynprof.sph.log 2>&1
done
