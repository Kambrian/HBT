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
	SUBCATALOGUE SubCat,SubCatB;
	SRCCATALOGUE SrcCat,SrcCatB;
	int Npro,*pro2dest;
	char subdirB[1024],srcdirB[1024];
	
	int Nsnap=86;
	int grpid=0,subid=29,i,j;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal",SUBCAT_DIR);	
	sprintf(subdirB,"%s/../patched_new",SUBCAT_DIR);
	sprintf(srcdirB,"%s/..",SUBCAT_DIR);
	logfile=stdout;
	if(argc!=3)
	{printf("usage: loadview [Nsnap], otherwise set Nsnap inside\n");fflush(stdout);}
	else
	{
	Nsnap=atoi(argv[1]);
	subid=atoi(argv[2]);
	}
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_src_catalogue(Nsnap,&SrcCat,SUBCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCatB,subdirB);
	load_src_catalogue(Nsnap,&SrcCatB,srcdirB);
	printf("A:\n%d,%d,%f\nCoM: %f,%f,%f\nVCoM:%f,%f,%f\n%f\n%f\n%f,%f,%f\n",SrcCat.SubLen[subid],SubCat.SubLen[subid],SrcCat.CoreFrac[subid],
		SubCat.Property[subid].CoM[0],SubCat.Property[subid].CoM[1],SubCat.Property[subid].CoM[2],
		SubCat.Property[subid].VCoM[0],SubCat.Property[subid].VCoM[1],SubCat.Property[subid].VCoM[2],
		SubCat.Property[subid].Pot,SubCat.Property[subid].Kin,
		SubCat.Property[subid].AM[0],SubCat.Property[subid].AM[1],SubCat.Property[subid].AM[2]);
	printf("B:\n%d,%d,%f\nCoM: %f,%f,%f\nVCoM:%f,%f,%f\n%f\n%f\n%f,%f,%f\n",SrcCatB.SubLen[subid],SubCatB.SubLen[subid],SrcCatB.CoreFrac[subid],
		SubCatB.Property[subid].CoM[0],SubCatB.Property[subid].CoM[1],SubCatB.Property[subid].CoM[2],
		SubCatB.Property[subid].VCoM[0],SubCatB.Property[subid].VCoM[1],SubCatB.Property[subid].VCoM[2],
		SubCatB.Property[subid].Pot,SubCatB.Property[subid].Kin,
		SubCatB.Property[subid].AM[0],SubCatB.Property[subid].AM[1],SubCatB.Property[subid].AM[2]);
		
	return 0;
}
