#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//to check for the binding energy of background particles

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	       
void match_sub(int IDArr[],int SubLen,int ID2Sub[],int Nsubs,int *desID,int *ShareLen);
int* prepare_id2sub(SUBCATALOGUE *A);
void free_id2sub(int *ID2H);

int main(int argc, char** argv)
{
	char catA_dir[1024]="/home/kambrain/data/8213/subcatS";
	char catB_dir[1024]="/home/kambrain/data/8213/subcat";
	char BName[1024]="HBT";
	char AName[1024]="SUBFIND";
	SUBCATALOGUE SubCatA,SubCatB;
	CATALOGUE Cat;
	int Nsnap=0;
	int i,grpid,subid,pid;
	int *ID2SubA,*ID2SubB;
	int *desID,*ShareLen;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=3)
	{printf("usage: %s [Nsnap][grpid]\n",argv[0]);exit(1);}
	
	Nsnap=atoi(argv[1]);
	grpid=atoi(argv[2]);
	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	#ifdef GRPINPUT_INDEX	
	for(i=0;i<Cat.Nids;i++)
		#ifdef PID_ORDERED
		Cat.PIDorIndex[i]++;//change from index to id
		#else
		Cat.PIDorIndex[i]=Pdat.PID[Cat.PIDorIndex[i]];
		#endif
	#endif
	load_subfind_catalogue(Nsnap,&SubCatA,catA_dir);
	load_sub_catalogue(Nsnap,&SubCatB,catB_dir);
	ID2SubA=prepare_id2sub(&SubCatA);
	ID2SubB=prepare_id2sub(&SubCatB);
	int ABack,AOther,ASub,ABackAll,*ALen;
	ABack=0;AOther=0;ASub=0;ABackAll=0;

	ALen=mymalloc(sizeof(int)*SubCatA.GrpLen_Sub[grpid]);
	for(i=0;i<SubCatA.GrpLen_Sub[grpid];i++)
	ALen[i]=0;
	for(i=Cat.Offset[grpid];i<Cat.Offset[grpid]+Cat.Len[grpid];i++)
	{
		pid=Cat.PIDorIndex[i];
		if(ID2SubB[pid]<0)//background in B
		{
			if(ID2SubA[pid]<0) //also background in A
				ABack++;
			else if(ID2SubA[pid]==SubCatA.GrpOffset_Sub[grpid])
				ASub++;
			else
			{
				AOther++;
				if(ID2SubA[pid]>=SubCatA.GrpOffset_Sub[grpid]+SubCatA.GrpLen_Sub[grpid])
				printf("%d,%d\n",pid,ID2SubA[pid]);
				ALen[ID2SubA[pid]-SubCatA.GrpOffset_Sub[grpid]]++;
			}
		}
		if(ID2SubA[pid]<0)
		ABackAll++;
	}
	ALen[0]=ASub;
	printf("back=%d,other=%d,central=%d,totalSFbk=%d\n",ABack,AOther,ASub,ABackAll);

	free_id2sub(ID2SubA);
	free_id2sub(ID2SubB);
	
	FILE *fp;
	char buf[1024],outdir[1024];
	sprintf(outdir,"%s/anal/MatchTo",catB_dir);
	mkdir(outdir,0755);
	sprintf(outdir,"%s/%s",outdir,AName);
	mkdir(outdir,0755);
	sprintf(buf,"%s/BackDistr_%d_%d",outdir,Nsnap,Nsnap);
	myfopen(fp,buf,"w");
	int sumlen=0;
	for(subid=0;subid<SubCatA.GrpLen_Sub[grpid];subid++)
	{
	fprintf(fp,"%d\t%d\t%d\n",subid,ALen[subid],SubCatA.SubLen[subid]);
	sumlen+=ALen[subid];
	}
	printf("%d\n",sumlen);
	fclose(fp);
	
	//~ erase_sub_catalogue(&SubCatA);
	//~ erase_sub_catalogue(&SubCatB);
	return 0;
}
#define CHECK_PID(pid)	if(pid<=NP_gas||pid>NP_gas+NP_DM) {printf("wrong pid:%d, out of range (%d,%d]\n",pid,NP_gas,NP_gas+NP_DM);	exit(1);}
int* prepare_id2sub(SUBCATALOGUE *A)
{
	int i,subid,pid;
	int *ID2H;
	
	ID2H=mymalloc(sizeof(int)*NP_DM);
	ID2H--;//move from id=0
	ID2H-=NP_gas;//move from gas
	#pragma omp parallel 
	{
	#pragma omp for
	for(i=1;i<=NP_DM;i++)//initialization
	{
		ID2H[i]=-1;/*return -1 if the PID does not belong to a Halo,
								i.e,we consider the backgroud as a halo with haloid=-1; 
								note that this only make sense when we try to find host for bound structures */
	}
	#pragma omp for private(i,pid,subid)
	for(subid=0;subid<A->Nsubs;subid++)
	{
		for(i=0;i<A->SubLen[subid];i++)
		{
			pid=A->PSubArr[subid][i];//PID ranges [1,NP_DM];
			CHECK_PID(pid)
			ID2H[pid]=subid;//haloIDs begins from id=0
		}
	}
	}
	return ID2H; 
}
void free_id2sub(int *ID2H)
{
	free(ID2H+1+NP_gas);
}
