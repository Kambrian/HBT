#!/bin/bash
export RUN_NUM=$1 #

export HBT_VERSION=8.7c
export OMP_SCHEDULE=dynamic
export OMP_DYNAMIC=FALSE # this works for icc v11.1.072, 
# but seems not always use the maximum available number of cpus allowed by NUM_THREADS
# so better not allow it to adjust. shit.
export KMP_BLOCKTIME=0  # this seems to help a lot, to reuse thread immediately after its job is finished
export OMP_NUM_THREADS=32
export OMP_NESTED=0

snapstart=0
snapend=59

date 

make clean
# make daily -e
make snap2 -e

echo =========================================================================

# ./HBT.$RUN_NUM $snapstart $snapend 

# ./haloprof.$RUN_NUM $snapend
# ./snaps.$RUN_NUM $snapend 0 $snapend
./snap2.$RUN_NUM $snapend 0 $snapend
# ./massfun_plot.$RUN_NUM $snapend 0  

# for((i=$(($snapend-1));i>=$snapstart;i-=1))
# do
#   echo "Snap=$i:" 
#   ./haloprof.$RUN_NUM $i 
# done
# chmod a-w ~/data/$RUN_NUM/subcat/profile/logbin/*

# ./NFW_fit.$RUN_NUM