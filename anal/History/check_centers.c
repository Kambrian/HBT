/*to check which is the best representation of subhalos' center*/
/*--Conclusion--*
% in general both the most-bound and the min-pot particle are good
% representation for the center, with min-pot particles have slightly larger dispersion;
% with only very rare occassions that
% min-pot particle can have substantial positional/velocity offset, i.e.,
% being ejected from the core.
% the velocity of min-pot particle also has low dispersion, ganranteed by
% the velocity dispersion profile (rho/sigma^3~r^alpha-->sigma~r^..., -->0
% as r-->0. i.e., the core is coherent and tight.
* ----*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

SUBCATALOGUE SubCat;
int minpot_id(int subid);

int main(int argc,char **argv)
{
	int i,subid,Nsnap,SnapNum;
	int *sp2pro;
	int Npro,Nsp;
	
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	
	if(argc!=3)
	{
		printf(" %s [Nsnap][subid]\n",argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);
	subid=atoi(argv[2]);
	
	FILE *fp,*fp2;
	char buf[1024];
	sprintf(outputdir,"%s/anal",SUBCAT_DIR);	
	mkdir(outputdir,0755);		
	logfile=stdout;
	
	sprintf(buf, "%s/check_centers_%d_%d",outputdir,Nsnap,subid);
	myfopen(fp,buf,"w");
	sprintf(buf, "%s/check_centers_%d_%d.final",outputdir,Nsnap,subid);
	myfopen(fp2,buf,"w");	
	#define PRNT3g(fd,p) fprintf(fd,",%g,%g,%g",p[0],p[1],p[2])
	
	int bid,pid,bidlist[MaxSnap]={},pidlist[MaxSnap]={},SnapS,SnapE;
	SnapS=Nsnap;
	while(subid>=0)
	{
		printf("%d,",Nsnap);fflush(stdout);
		load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
		load_particle_data(Nsnap,SNAPSHOT_DIR);
		fresh_ID2Index(&SubCat,FRSH_SUBCAT);
		
		bid=SubCat.PSubArr[subid][0];
		pid=minpot_id(subid);
		bidlist[Nsnap]=Pdat.PID[bid];
		pidlist[Nsnap]=Pdat.PID[pid];
		fprintf(fp,"%d,%d,%d,%d,%d",Nsnap,subid,SubCat.SubRank[subid],Pdat.PID[bid],Pdat.PID[pid]);
		PRNT3g(fp,SubCat.Property[subid].CoM);
		PRNT3g(fp,SubCat.Property[subid].VCoM);
		PRNT3g(fp,Pdat.Pos[bid]);
		PRNT3g(fp,Pdat.Pos[pid]);
		PRNT3g(fp,Pdat.Vel[bid]);
		PRNT3g(fp,Pdat.Vel[pid]);
		fprintf(fp,"\n");
		
		subid=SubCat.HaloChains[subid].ProSubID;
		load_sp2pro(Nsnap,&Npro,&Nsp,&sp2pro,SUBCAT_DIR);
		if(subid>=Npro)
		{
			printf("splinter: %d,%d\n",Nsnap-1,subid);
			subid=sp2pro[subid];
		 }
		free_sp2pro(sp2pro,Npro,Nsp);
		free_sub_catalogue(&SubCat);
		Nsnap--;
	}
	fclose(fp);
	puts("\n");
	/*==pos of the central particles at final time===*/
	load_particle_data(SnapS,SNAPSHOT_DIR);
	fresh_ID2Index(bidlist,MaxSnap);
	fresh_ID2Index(pidlist,MaxSnap);
	SnapE=Nsnap;
	for(Nsnap=SnapS;Nsnap>SnapE;Nsnap--)
	{
	fprintf(fp2,"%d",Nsnap);
	bid=bidlist[Nsnap];
	pid=pidlist[Nsnap];
	PRNT3g(fp2,Pdat.Pos[bid]);
	PRNT3g(fp2,Pdat.Pos[pid]);
	PRNT3g(fp2,Pdat.Vel[bid]);
	PRNT3g(fp2,Pdat.Vel[pid]);
	fprintf(fp2,"\n");
	}
	fclose(fp2);
	return 0;
}

int minpot_id(int subid)
{
//find particle index of min-pot 	
int ipotmin,i,pid;
double potmin,*pot;

	potmin=0.;ipotmin=0;
	tree_tree_allocate(TREE_ALLOC_FACTOR*SubCat.SubLen[subid],SubCat.SubLen[subid]);
	maketree(SubCat.SubLen[subid],SubCat.PSubArr[subid],Pdat.Pos);
	pot=mymalloc(sizeof(double)*SubCat.SubLen[subid]);
	for(i=0;i<SubCat.SubLen[subid];i++)
	{
		pot[i]=tree_treeevaluate_potential(Pdat.Pos[SubCat.PSubArr[subid][i]],SubCat.PSubArr[subid],Pdat.Pos);
		if(potmin>=pot[i])
		{
		ipotmin=i;
		potmin=pot[i];
		}
	}
	//~ printf("%d:\t%d\t%f\n",subid,ipotmin,potmin);
	pid=SubCat.PSubArr[subid][ipotmin];
	free(pot);
	tree_tree_free();
	return pid;
}
