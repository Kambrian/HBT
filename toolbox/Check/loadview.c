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
	HBTInt IDArr[10];
	
	HBTInt Nsnap=0;
	HBTInt grpid,subid,pid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],(int)Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
// 	load_sub_table(Nsnap, &SubCat, SUBCAT_DIR);
 load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
 load_particle_header(Nsnap, SNAPSHOT_DIR);
	//~ load_src_catalogue(Nsnap,&SrcCat,SUBCAT_DIR);
 /*
	if(Nsnap<MaxSnap-1)
 load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);
 load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);*/
	load_particle_data_bypart(Nsnap,SNAPSHOT_DIR, FLAG_LOAD_ID);
	for(pid=0;pid<10;pid++) IDArr[pid]=Pdat.PID[pid];
	fill_PIDHash();
	for(pid=0;pid<10;pid++)
	{
	  printf("%ld, %ld, %ld, %ld\n", pid, lookup_ID2Ind(Pdat.PID[pid]), Pdat.PID[pid], Pdat.PID[lookup_ID2Ind(Pdat.PID[pid])]);
	}
	fresh_ID2Index(IDArr, 10);
	fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);	//fresh_ID2Index(&SrcCat,-3);
	free_PIDHash();
	free_particle_data();
	load_particle_data_bypart(Nsnap, SNAPSHOT_DIR, FLAG_LOAD_POS);
	
	HBTInt i,j=0,k=0;
	for(i=0;i<SubCat.Ngroups;i++)
	{
		subid=SubCat.GrpOffset_Sub[i];
		if(0==SubCat.SubLen[subid]) j++;
		if(SubCat.SubLen[subid]>=exp10(0.5)/header.mass[1]) k++;
	}
	printf(""HBTIFMT" halos, "HBTIFMT" subhalos loaded;"
	       " "HBTIFMT" halos are not bound;"
		   " "HBTIFMT" central subhalos have M>1e10.5 Msun/h\n",
		   SubCat.Ngroups,SubCat.Nsubs,j,k);
	
	printf("Halo #0 Mass: "HBTIFMT"\n",Cat.Len[0]);
	printf("Halo #0 has "HBTIFMT" subhalos\n", SubCat.GrpLen_Sub[0]);
	printf("Subhalo #0 Mass: "HBTIFMT"\n",SubCat.SubLen[0]);
	printf("Subhalo #0 Host halo ID: "HBTIFMT"\n", SubCat.HaloChains[0].HostID);
	printf("Subhalo #0's progenitor subhalo ID:"HBTIFMT"\n", SubCat.HaloChains[0].ProSubID);
	subid=SubCat.GrpOffset_Sub[1];
	printf("Halo #1's central subhalo's center: %f, %f, %f\n", SubCat.Property[subid].CoM[0],SubCat.Property[subid].CoM[1],SubCat.Property[subid].CoM[2]);
// 	pid=SubCat.PSubArr[subid][0];
// 	printf("Halo #1's central subhalo's most-bound position: %f, %f, %f\n", Pdat.Pos[pid][0],Pdat.Pos[pid][1],Pdat.Pos[pid][2]);
	
	free_catalogue(&Cat);
	erase_sub_catalogue(&SubCat);
// 	erase_src_catalogue(&SrcCat);
// 	if(Nsnap<MaxSnap-1)
// 	free_pro2dest(pro2dest);
// 	free_sp2pro(sp2pro,Npro,Nsplitter);
	return 0;
}
