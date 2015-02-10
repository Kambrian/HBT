#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc, char** argv)
{
	SUBCATALOGUE SubCat;
	
	HBTInt Nsnap=MaxSnap-1;
	HBTInt subid=0;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=3)
	{printf("usage: %s [Nsnap] [subid], otherwise Nsnap=%d, subid=%d\n",argv[0],(int)Nsnap, (int)subid);fflush(stdout);}
	else
	{
	  Nsnap=atoi(argv[1]);
	  subid=atoi(argv[2]);
	}
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&SubCat,-2);	//fresh_ID2Index(&SrcCat,-3);
	free_PIDHash();

    #define READPOS(x) x[0],x[1],x[2]
	printf("Cen=%g,%g,%g\n", READPOS(SubCat.Property[subid].CoM));
	printf("VCen=%g,%g,%g\n", READPOS(SubCat.Property[subid].VCoM));
	printf("MostBound Cen: %g,%g,%g\n",READPOS(Pdat.Pos[SubCat.PSubArr[subid][0]]));
	printf("MostBound VCen: %g,%g,%g\n",READPOS(Pdat.Vel[SubCat.PSubArr[subid][0]]));
	
	erase_sub_catalogue(&SubCat);
	return 0;
}
