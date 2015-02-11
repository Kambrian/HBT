#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//to match subhalos from catB to catA;

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
	char catB_dir[1024]="/home/kambrain/data/6702DM/subcat";
	char catA_dir[1024]="/home/kambrain/data/6702DM/subcatS";
	char AName[1024]="SUBFIND";
	char BName[1024]="HBT";
	SUBCATALOGUE SubCatA,SubCatB;
	int NsnapA,NsnapB;
	int grpid,subid,pid;
	int *ID2SubA,*ID2SubB;
	int *desID,*ShareLen;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=3)
	{printf("usage: %s [NsnapA],[NsnapB]\n",argv[0]);fflush(stdout);exit(1);}
	else
	{
	NsnapA=atoi(argv[1]);
	NsnapB=atoi(argv[2]);
	}

	load_sub_catalogue(NsnapA,&SubCatA,catA_dir);
	load_sub_catalogue(NsnapB,&SubCatB,catB_dir);
	ID2SubA=prepare_id2sub(&SubCatA);
	//~ ID2SubB=prepare_id2sub(&SubCatB);
	desID=mymalloc(sizeof(int)*SubCatB.Nsubs);
	ShareLen=mymalloc(sizeof(int)*SubCatB.Nsubs);
	#pragma omp parallel for
	for(subid=0;subid<SubCatB.Nsubs;subid++)
	match_sub(SubCatB.PSubArr[subid],SubCatB.SubLen[subid],ID2SubA,SubCatA.Nsubs,desID+subid,ShareLen+subid);
	
	free_id2sub(ID2SubA);
	
	FILE *fp;
	char buf[1024],outdir[1024];
	sprintf(outdir,"%s/anal/MatchTo",catB_dir);
	mkdir(outdir,0755);
	sprintf(outdir,"%s/%s",outdir,AName);
	mkdir(outdir,0755);
	sprintf(buf,"%s/%d_%d",outdir,NsnapB,NsnapA);
	myfopen(fp,buf,"w");
	for(subid=0;subid<SubCatB.Nsubs;subid++)
	fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",subid,desID[subid],ShareLen[subid],SubCatB.SubLen[subid],desID[subid]<0?-1:SubCatA.SubLen[desID[subid]]);
	fclose(fp);
	
	//~ erase_sub_catalogue(&SubCatA);
	//~ erase_sub_catalogue(&SubCatB);
	return 0;
}
void match_sub(int IDArr[],int SubLen,int ID2Sub[],int Nsubs,int *desID,int *ShareLen)
{
	int i,subid;
	int *LinkCount;
	float *LinkWght;
	LinkCount=mymalloc(sizeof(int)*(Nsubs+1)); //include id=-1, the background
	LinkWght=mymalloc(sizeof(float)*(Nsubs+1));
	LinkCount++;
	LinkWght++;
	for(subid=-1;subid<Nsubs;subid++)
	{
		 LinkCount[subid]=0;
		 LinkWght[subid]=0.;
	 }
	for(i=0;i<SubLen;i++)
	{
		subid=ID2Sub[IDArr[i]];
		LinkCount[subid]++;
		LinkWght[subid]+=powf(i+1,-0.6);
	}
	subid=Fmax_of_vec(LinkWght-1,Nsubs+1)-1;
	//subid=Fmax_of_vec(LinkWght,Nsubs);//select only within real subs
	//~ if(LinkCount[subid])
	*desID=subid;
	//~ else
	//~ *desID=-1;
	*ShareLen=LinkCount[subid];
	myfree(LinkCount-1);
	myfree(LinkWght-1);
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