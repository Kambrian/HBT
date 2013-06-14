#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "intra_vars.h"
#include "iovars.h"
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"

int main(int argc,char **argv)
{
char gasdir[1024]=GASCAT_DIR;
char buf[1024],filemode[4];
int Nsnap,i,SnapRange[2];
GASHALOCAT GCat;
GASSUBCAT GSubCat;
GASSRCCAT GSrcCat;
SUBCATALOGUE SubCat;

while(1)
{
printf("Nsnap=");
scanf("%d",&Nsnap);getchar();
load_gashalocat(Nsnap,&GCat,GASCAT_DIR);
load_gassubcat(Nsnap,&GSubCat,GASCAT_DIR);
load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
printf("cat loaded,press any key to continue\n");
getchar();
erase_gassubcat(&GSubCat);
erase_sub_catalogue(&SubCat);
}

return 0;
}
