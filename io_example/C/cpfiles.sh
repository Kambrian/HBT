#!/bin/bash
cp ../../anal/Check/loadview.c .
cd library
cp ../../../gadget_io -R .
cp ../../../JING_io -R .
cp ../../../intra_vars.c  .
cp ../../../intra_vars.h  .
cp ../../../datatypes.h  .
cp ../../../mymath.c  .
cp ../../../proto.h  .
cp ../../../sub_IO.c  .
cp ../../../Makefile.runs  .
cp ../../../param -R .
cd history
cp ../../../../history/history_io.c 
cp ../../../../history/history_vars.h 
cp ../../../../history/history_proto.h 

