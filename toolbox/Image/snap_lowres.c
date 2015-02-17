#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "param6702low.h"
#include "allvars.h"
#include "proto.h"

#define NGRID 1500
int mapxy[NGRID][NGRID]={0},mapxz[NGRID][NGRID]={0},mapyz[NGRID][NGRID]={0};

int main()
{
	int Nsnap,i,j,pid;
	FILE *fp;
	char buf[1024];
	float step[3],range[3][2]={{147873.375000,154674.093750},{168185.203125,175743.796875},{144266.937500,149686.625000}};
	int grid[3];
	char outputdir[1024]="/home/kambrain/bt6/anal/low6702";
	
	logfile=stdout;
	Nsnap=99;

	load_particle_data(Nsnap,SNAPSHOT_DIR);
	for(i=0;i<3;i++)
		step[i]=(range[i][1]-range[i][0])/NGRID;
	printf("size: %f,%f,%f\n",(range[0][1]-range[0][0]),(range[1][1]-range[1][0]),(range[2][1]-range[2][0]));
 //==================2d map======================//
	for(i=0;i<NP_DM;i++)
	{
		pid=i;
		for(j=0;j<3;j++)
		{
			if(Pdat.Pos[pid][j]>=range[j][0]&&Pdat.Pos[pid][j]<=range[j][1])
			{
			grid[j]=floor((Pdat.Pos[pid][j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
			}
		}
		mapxy[grid[0]][grid[1]]++;
		mapxz[grid[0]][grid[2]]++;
		mapyz[grid[1]][grid[2]]++;
	}
	
	sprintf(buf,"%s/fofmapxy_%03d",outputdir,Nsnap);
	if((fp=fopen(buf,"w"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		exit(1);
	}
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxy[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/fofmapxz_%03d",outputdir,Nsnap);
	if((fp=fopen(buf,"w"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		exit(1);
	}
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/fofmapyz_%03d",outputdir,Nsnap);
	if((fp=fopen(buf,"w"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		exit(1);
	}
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapyz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
				
return 0;
}
		
