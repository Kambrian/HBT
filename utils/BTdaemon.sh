#!/bin/bash
while true
do
if pgrep BT.6113
then
	sleep 30m
else
	./BTrun
	echo BT running from `date`
	exit 0
fi
done
