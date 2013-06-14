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

#include "gas_stuff.h"


int main()
{
char outputdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5/GasCat/nonthermal_exclusive";
char inputdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5";
char fofdir[512]="/SANdisk5/kambrain/Sim6702/FoFCat"; //"/home/kambrain/fof_hy";
char snapdir[512]="/SANdisk4/data/NewR/SIM6702";
	
int Nsnap=99;
CATALOGUE Cat;
SUBCATALOGUE SubCat;
GASCATALOGUE GasCat;
int haloid,subid,suboffset;
char buf[1024];

	//logfile=stdout;	
	sprintf(buf,"%s/logfile_gas",outputdir);
	if(!(logfile=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}	
	
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	load_particle_data(Nsnap,&Cat,&SubCat,snapdir);
	load_gas_data(Nsnap,snapdir);
	create_gas_cat(SubCat.Nsubs,&GasCat);
	Rtidal=mymalloc(sizeof(float)*SubCat.Nsubs);
	load_tidal_radius(Nsnap,Rtidal,SubCat.Nsubs,inputdir);
	
	makell(Gdat.Pos,NP_gas);
	
	for(haloid=0;haloid<Cat.Ngroups;haloid++)
	{
		unbind_gas_recursive(SubCat.GrpOffset_Sub[haloid],&SubCat,&GasCat);
	}
	for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++)
	{
	collect_gas_particles(subid,&SubCat,&GasCat);
	unbindgas(GasCat.SubLen+subid,GasCat.PSubArr+subid,SubCat.SubLen[subid],SubCat.PSubArr[subid]);
	}
	
	suboffset=0;
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		GasCat.SubOffset[subid]=suboffset;
		suboffset+=GasCat.SubLen[subid];
	}
	GasCat.Nids=suboffset;
	
	save_gas_cat(Nsnap,&GasCat,outputdir);
	
	return 0;
}
 
