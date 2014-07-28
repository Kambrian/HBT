#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define NO_SUBCAT
// #define SUBFIND_DIR "/home/jvbq85/data/HBT/data/AqE5W/subfind"
// #define OUTDIR "/home/jvbq85/data/HBT/data/AqE5W/subcat/anal/subfind"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
//#define load_group_catalogue load_group_catalogue_v3
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

#define NGRID 500
	int mapxy[NGRID][NGRID]={0},mapxz[NGRID][NGRID]={0},mapyz[NGRID][NGRID]={0};

int main(int argc,char **argv)
{
	CATALOGUE Cat;
#ifndef NO_SUBCAT
	SUBCATALOGUE SubCat;
#endif
	int SnapLoad,SnapPlot,i,j,pid,subid,grpid;
	FILE *fp,*fpr,*fpv;
	char buf[1024];
	float step[3],range[3][2];

	int *PIndex,grid[3],ipotmin=0;
	double *pot,potmin=0.,com[3];
	char outputdir[1024];
	#ifdef SUBFIND_DIR
	sprintf(outputdir,"%s/image",OUTDIR);	
	#else
	sprintf(outputdir,"%s/anal/image",SUBCAT_DIR);
	#endif
	mkdir(outputdir,0755);
	logfile=stdout;

	if(argc!=4)
	{
		printf("usage: %s [SnapLoad] [grpid] [SnapPlot], otherwise set Nsnap inside\n",argv[0]);
		fflush(stdout);
	}
	else
	{
		SnapLoad=atoi(argv[1]);
		grpid=atoi(argv[2]);
		SnapPlot=atoi(argv[3]);
	}

	load_group_catalogue(SnapLoad,&Cat,GRPCAT_DIR);
#ifndef NO_SUBCAT
	load_sub_catalogue(SnapLoad,&SubCat,SUBCAT_DIR);
#endif
	load_particle_data(SnapPlot,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&Cat,-1); 	
#ifndef NO_SUBCAT
	fresh_ID2Index(&SubCat,-2);	
#endif
	free_PIDHash();
	
	PIndex=Cat.PIDorIndex+Cat.Offset[grpid];
	for(i=0;i<3;i++)
		for(j=0;j<2;j++)
			range[i][j]=Pdat.Pos[PIndex[0]][i];
	for(i=1;i<Cat.Len[grpid];i++)
	{
		pid=PIndex[i];
		for(j=0;j<3;j++)
		{
		if(Pdat.Pos[pid][j]<range[j][0])
			range[j][0]=Pdat.Pos[pid][j];
		else if(Pdat.Pos[pid][j]>range[j][1])
			range[j][1]=Pdat.Pos[pid][j];
		}
	}
	for(i=0;i<3;i++)
		step[i]=(range[i][1]-range[i][0])/NGRID;
	sprintf(buf,"%s/fofsize_%03d_%d.%03d",outputdir,SnapLoad,grpid,SnapPlot);
	myfopen(fp,buf,"w");
	for(i=0;i<3;i++)	
		fprintf(fp,"%f,%f\n",range[i][0],range[i][1]);
	fclose(fp);
	//~ //==================2d map======================//
	for(i=0;i<Cat.Len[grpid];i++)
	{
		pid=PIndex[i];
		for(j=0;j<3;j++)
		{
			grid[j]=floor((Pdat.Pos[pid][j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
			if(grid[j]<0) grid[j]=0;
		}
		mapxy[grid[0]][grid[1]]++;
		mapxz[grid[0]][grid[2]]++;
		mapyz[grid[1]][grid[2]]++;
	}
	
	sprintf(buf,"%s/fofmapxy_%03d_%d.%03d",outputdir,SnapLoad,grpid,SnapPlot);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxy[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/fofmapxz_%03d_%d.%03d",outputdir,SnapLoad,grpid,SnapPlot);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/fofmapyz_%03d_%d.%03d",outputdir,SnapLoad,grpid,SnapPlot);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapyz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	
#ifndef NO_SUBCAT	
	//===========most bound===============//
	sprintf(buf,"%s/subcenmap_%03d_%d.%03d",outputdir,SnapLoad,grpid,SnapPlot);
	myfopen(fp,buf,"w");
	for(i=0;i<SubCat.GrpLen_Sub[grpid];i++)
	{
		subid=SubCat.GrpOffset_Sub[grpid]+i;
		pid=SubCat.PSubArr[subid][0];
		for(j=0;j<3;j++)
		{
			grid[j]=floor((Pdat.Pos[pid][j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
			fprintf(fp,"%d\t",grid[j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);	
	sprintf(buf,"%s/subcen_%03d_%d.%03d",outputdir,SnapLoad,grpid,SnapPlot);
	myfopen(fpr,buf,"w");
		for(i=0;i<SubCat.GrpLen_Sub[grpid];i++)
	{
		subid=SubCat.GrpOffset_Sub[grpid]+i;
		pid=SubCat.PSubArr[subid][0];
		fprintf(fpr,"%f\t%f\t%f\n",Pdat.Pos[pid][0],Pdat.Pos[pid][1],Pdat.Pos[pid][2]);
	}
	fclose(fpr);
	
	//=================min of pot==================//
	//~ sprintf(buf,"%s/subminmap_%03d_%d",outputdir,Nsnap,grpid);
	//~ if((fp=fopen(buf,"w"))==NULL)
	//~ {
		//~ fprintf(logfile,"error: file open failed for %s!\n",buf);
		//~ exit(1);
	//~ }
	//~ sprintf(buf,"%s/submin_%03d_%d",outputdir,Nsnap,grpid);
	//~ if((fpr=fopen(buf,"w"))==NULL)
	//~ {
		//~ fprintf(logfile,"error: file open failed for %s!\n",buf);
		//~ exit(1);
	//~ }
	//~ for(subid=0;subid<SubCat.GrpLen_Sub[grpid];subid++)
	//~ {
		//~ potmin=0.;ipotmin=0;
		//~ tree_tree_allocate(TREE_ALLOC_FACTOR*SubCat.SubLen[subid],SubCat.SubLen[subid]);
		//~ maketree(SubCat.SubLen[subid],SubCat.PSubArr[subid]);
		 //~ pot=mymalloc(sizeof(double)*SubCat.SubLen[subid]);
		//~ for(i=0;i<SubCat.SubLen[subid];i++)
		//~ {
			//~ pot[i]=tree_treeevaluate_potential(SubCat.PSubArr[subid][i],SubCat.PSubArr[subid]);
			//~ if(potmin>=pot[i])
			//~ {
			//~ ipotmin=i;
			//~ potmin=pot[i];
			//~ }
		//~ }
		//~ printf("%d:\t%d\t%f\n",subid,ipotmin,potmin);
		//~ pid=SubCat.PSubArr[subid][ipotmin];
		//~ for(j=0;j<3;j++)
		//~ {
			//~ grid[j]=floor((Pdat.Pos[pid][j]-range[j][0])/step[j]);
			//~ if(grid[j]>=NGRID) grid[j]=NGRID-1;
			//~ fprintf(fp,"%d\t",grid[j]);
		//~ }
		//~ fprintf(fp,"\n");
		//~ fprintf(fpr,"%f\t%f\t%f\n",Pdat.Pos[pid][0],Pdat.Pos[pid][1],Pdat.Pos[pid][2]);
		//~ free(pot);
		//~ tree_tree_free();
	//~ }
	//~ fclose(fp);
	//~ fclose(fpr);
	
	//========================com=================//
	sprintf(buf,"%s/subcommap_%03d_%d",outputdir,SnapLoad,grpid);
	myfopen(fp,buf,"w");
	sprintf(buf,"%s/subcom_%03d_%d",outputdir,SnapLoad,grpid);
	myfopen(fpr,buf,"w");
	for(i=0;i<SubCat.GrpLen_Sub[grpid];i++)//NOTE: if not first fof, some modification for subid
	{
		subid=SubCat.GrpOffset_Sub[grpid]+i;
		for(j=0;j<3;j++)
		{	
			com[j]=SubCat.Property[subid].CoM[j];
			grid[j]=floor((com[j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
			fprintf(fp,"%d\t",grid[j]);
		}
		fprintf(fp,"\n");
		fprintf(fpr,"%d\t%f\t%f\t%f\n",SubCat.SubLen[subid],com[0],com[1],com[2]);
	}
	fclose(fp);
	fclose(fpr);
#endif
				
return 0;
}
