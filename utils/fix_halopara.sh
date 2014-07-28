#!/bin/sh
#usage: fix_halopara.sh -i param*.h
sed '
/\#define HALO_PARA/ c\
\#ifndef DISABLE_HALO_PARA\
\t\#define HALO_PARA\
\#endif
' $@