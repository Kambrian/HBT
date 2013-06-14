#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"

int main(int argc,char **argv)
{
	char buf[1024];
	CATALOGUE Cat;
	GASHALOCAT GCat;
	HBTInt i,pid,hostid,Nsnap,*DMPIndex,Offset;
	HBTReal hguess;
	time_t time_p[10];
	
	sprintf(buf,"%s/%s.gashalo",GASCAT_DIR,LOGFILE_NAME);
	if((logfile=fopen(buf,"w"))==NULL)	{fprintf(stderr,"Error opening file %s\n",buf);exit(1);}
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
	time_p[0]=time(NULL);
	load_particle_data(Nsnap,SNAPSHOT_DIR);	
	load_gas_data(Nsnap,SNAPSHOT_DIR);
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	#ifndef GRPINPUT_INDEX
	fresh_ID2Index(&Cat,-1); 	//fofcat of JING's data are originally PIndex rather than PID
	#endif
	prepare_ind2halo(&Cat);
	
	fprintf(logfile,"making tree ...\n");fflush(logfile);time_p[1]=time(NULL);
	tree_tree_allocate(TREE_ALLOC_FACTOR*NP_DM,NP_DM);
	DMPIndex=mymalloc(sizeof(int)*NP_DM);
	for(i=0;i<NP_DM;i++)
		DMPIndex[i]=i;
	maketree(NP_DM,DMPIndex,Pdat.Pos);
		
	fprintf(logfile,"finding nearest dm neighbor ...\n");fflush(logfile);time_p[2]=time(NULL);
	GCat.ID2Halo=mymalloc(sizeof(int)*NP_GAS);
	hguess=guess_ngb_range(1);
	#pragma omp parallel for private(i,pid) schedule(dynamic)
	for(i=0;i<NP_GAS;i++)
	{
		pid=treesearch_nearest(Gdat.Pos[i],hguess,DMPIndex,Pdat.Pos);
		GCat.ID2Halo[i]=Cat.ID2Halo[pid];
		//~ printf("%d\n",i);
	}
	tree_tree_free();
	free(DMPIndex);
	free_catalogue(&Cat);
		
	fprintf(logfile,"building halo access info...\n");fflush(logfile);time_p[3]=time(NULL);
	GCat.Ngroups=Cat.Ngroups;
	GCat.Len=mymalloc(sizeof(int)*GCat.Ngroups);
	GCat.Offset=mymalloc(sizeof(int)*GCat.Ngroups);
	for(i=0;i<GCat.Ngroups;i++)
		GCat.Len[i]=0;
	for(i=0;i<NP_GAS;i++)
		if(GCat.ID2Halo[i]>=0)
			GCat.Len[GCat.ID2Halo[i]]++;	
	Offset=0;
	for(i=0;i<GCat.Ngroups;i++)
	{
		GCat.Offset[i]=Offset;
		Offset+=GCat.Len[i];
	}
	GCat.Nids=Offset;
	
	fprintf(logfile,"gathering halo particles ...\n");fflush(logfile);time_p[4]=time(NULL);
	GCat.PIDorIndex=mymalloc(sizeof(int)*GCat.Nids);
	for(i=0;i<GCat.Ngroups;i++)
		GCat.Len[i]=0;
	for(i=0;i<NP_GAS;i++)
	{
		if((hostid=GCat.ID2Halo[i])>=0)
		{			
			GCat.PIDorIndex[GCat.Offset[hostid]+GCat.Len[hostid]]=Gdat.PID[i];
			GCat.Len[hostid]++;
		}
	}
	
	fprintf(logfile,"saving halo ...\n");fflush(logfile);time_p[5]=time(NULL);
	save_gashalocat(Nsnap,&GCat,GASCAT_DIR);
	free_gashalocat(&GCat);
	free(GCat.ID2Halo);
	free_gas_data();
	time_p[6]=time(NULL);
	fprintf(logfile,"program timeing:\nload:%ld sec\ttreebuild: %ld sec\tngbsearch: %ld sec\thalostat: %ld sec\thalofill: %ld sec\tsave: %ld sec\n",
		(time_p[1]-time_p[0]),(time_p[2]-time_p[1]),(time_p[3]-time_p[2]),(time_p[4]-time_p[3]),(time_p[5]-time_p[4]),(time_p[6]-time_p[5]));
	fprintf(logfile,"total: %ld min\n",(time_p[6]-time_p[0])/60);
	}
	return 0;
}
