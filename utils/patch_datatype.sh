#!/bin/bash
dir=$1

cd $dir

for file in `ls *.c`
do
i=`grep '\#include "datatypes.h"' $file -c`

if [ $i -eq 0 ]; then
echo fixing $file
addtype.sh <$file >$file.fix
mv $file.fix $file
#equivalent,in-place edit with -i: addtype.sh -i $file
else
echo $file ok
fi
done
