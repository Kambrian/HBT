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
	CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	SRCCATALOGUE SrcCat;
	HBTInt Nsubs,*pro2dest,Npro,Nsplitter,*sp2pro;
	
	HBTInt Nsnap=0;
	HBTInt grpid,subid,pid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_src_catalogue(Nsnap,&SrcCat,SUBCAT_DIR);
	if(Nsnap<MaxSnap-1)
	load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);
	load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);	fresh_ID2Index(&SrcCat,-3);
	free_PIDHash();
	HBTInt i,j=0,k=0;
	for(i=0;i<SubCat.Ngroups;i++)
	{
		subid=SubCat.GrpOffset_Sub[i];
		if(0==SubCat.SubLen[subid]) j++;
		if(SubCat.SubLen[subid]>=exp10(0.5)/header.mass[1]) k++;
	}
	printf("%d halos, %d subhalos loaded;"
	       " %d halos are not bound;"
		   " %d central subhalos have M>1e10.5 Msun/h\n",
		   SubCat.Ngroups,SubCat.Nsubs,j,k);
	
	printf("Halo #0 Mass: %d\n",Cat.Len[0]);
	printf("Halo #0 has %d subhalos\n", SubCat.GrpLen_Sub[0]);
	printf("Subhalo #0 Mass: %d\n",SubCat.SubLen[0]);
	printf("Subhalo #0 Host halo ID: %d\n", SubCat.HaloChains[0].HostID);
	printf("Subhalo #0's progenitor subhalo ID:%d\n", SubCat.HaloChains[0].ProSubID);
	subid=SubCat.GrpOffset_Sub[1];
	printf("Halo #1's central subhalo's center: %f, %f, %f\n", SubCat.Property[subid].CoM[0],SubCat.Property[subid].CoM[1],SubCat.Property[subid].CoM[2]);
	pid=SubCat.PSubArr[subid][0];
	printf("Halo #1's central subhalo's most-bound position: %f, %f, %f\n", Pdat.Pos[pid][0],Pdat.Pos[pid][1],Pdat.Pos[pid][2]);
	
	free_catalogue(&Cat);
	erase_sub_catalogue(&SubCat);
	erase_src_catalogue(&SrcCat);
	if(Nsnap<MaxSnap-1)
	free_pro2dest(pro2dest);
	free_sp2pro(sp2pro,Npro,Nsplitter);
	return 0;
}
