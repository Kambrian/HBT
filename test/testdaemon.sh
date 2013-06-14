#!/bin/bash
while true
do
if pgrep BT.6113
then
	echo BT.6113 is there with pid=`pgrep BT.6113`
	sleep 5s
#	exit 1
else
	echo BT not running from `date`
	exit 0
fi
done
