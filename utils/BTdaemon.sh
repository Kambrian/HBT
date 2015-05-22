#!/bin/bash
while true
do
if pgrep BT.6113
then
	sleep 30m
else
# 	./BTrun
# 	echo BT running from `date`
	echo at `date` |mail -s "HBT finished" hanjiaxin@gmail.com -- -f jxhan@uv2000.shao.ac.cn
	exit 0
fi
done
