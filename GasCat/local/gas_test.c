#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "../globals.c"
#include "../mymath.c"
#include "../load_group.c"
#include "../tree.c"

#include "gas_stuff.c"


int main()
{
char outputdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5";
char inputdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5/GasCat/";
char fofdir[512]="/SANdisk5/kambrain/Sim6702/FoFCat"; //"/home/kambrain/fof_hy";
char snapdir[512]="/SANdisk4/data/NewR/SIM6702";
char gaskind[512];
	
int Nsnap=99;
//~ CATALOGUE Cat;
//~ SUBCATALOGUE SubCat;
GASCATALOGUE GasCat;
int haloid,subid,suboffset;
char buf[1024];
int i;

logfile=stdout;	
printf("input gas subdir:\n");
scanf("%s,%d\n",gaskind,&i);
strcat(inputdir,gaskind);

	load_gas_cat(Nsnap,&GasCat,inputdir);
	for(subid=0;subid<10;subid++)
	printf("%d\n",GasCat.SubLen[subid]);
	
	return 0;
}
 
