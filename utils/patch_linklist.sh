#!/bin/bash
dir=$1

cd $dir

for file in `ls *.c`
do
i=`grep '\#include "linkedlist.*.c"' $file -c`

if [ $i -eq 0 ]; then
echo $file ok
else
echo fixing $file
sed 's/\#include "\(linkedlist.*\).c"/\#include "\1.h"/' <$file >$file.fix
#addtype.sh <$file >$file.fix
mv $file.fix $file
fi
done
