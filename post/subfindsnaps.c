#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "globals.c"
#include "mymath.c"
#include "load_group.c"
#include "tree.c"
#include "subfindread.c"
//#include "binding.c"

#define NGRID 1500
	//~ int mapxy[NGRID][NGRID]={0},mapxz[NGRID][NGRID]={0},mapyz[NGRID][NGRID]={0};
	
int main()
{
	CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	int Nsnap,i,j,pid,subid,grpid;
	FILE *fp,*fpr,*fpv;
	char buf[1024];
	float step[3],range[3][2];

	int *PIndex,*PIndex_all,grid[3],ipotmin=0;
	double *pot,potmin=0.,com[3];
	char inputdir[512]="/SANdisk1/wenting/mergetree/SIM6702/subcatS";
	char fofdir[512]="/raid1/hywang/ReSim/SIM6702/group_catalogue"; //"/home/kambrain/fof_hy";
	char snapdir[512]="/SANdisk4/data/NewR/SIM6702";
	char outputdir[1024]="/SANdisk5/kambrain/Sim6702/SubCat5/anal";
	
	logfile=stdout;
	Nsnap=99;
	grpid=0;
	subid=3;

	load_group_catalogue(Nsnap,&Cat,fofdir);
	subread(Nsnap,inputdir);
	subfindcat2BT(&SubCat);
	load_particle_data(Nsnap,&Cat,&SubCat,snapdir);
		
	PIndex_all=mymalloc(sizeof(int)*NP_SIM);
	PIndex_all=PIndex_all-1;/*shift PIndex by -1 element so that 
					it is accessed through PIndex[1~N],
					i.e.,its subscript ranges the same as particle ID range.*/
	for(i=1;i<=NP_SIM;i++)
		PIndex_all[i]=-1;
	for(i=0;i<NP_DM;i++)
		PIndex_all[Pdat.PID[i]]=i;//suppose PID begins from ID=1;
		
	//~ sprintf(buf,"%s/subpos_%03d_%03d",outputdir,Nsnap,subid);
	//~ if((fpr=fopen(buf,"w"))==NULL)
	//~ {
		//~ fprintf(logfile,"error: file open failed for %s!\n",buf);
		//~ exit(1);
	//~ }
	
	//~ for(i=0;i<SubCat.SubLen[subid];i++)
	//~ {
		//~ pid=SubCat.PSubArr[subid][i];
		//~ fprintf(fpr,"%f\t%f\t%f\n",Pdat.Pos[pid][0],Pdat.Pos[pid][1],Pdat.Pos[pid][2]);
	//~ }
	//~ fclose(fpr);
			
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
	printf("size: %f,%f,%f\n",(range[0][1]-range[0][0]),(range[1][1]-range[1][0]),(range[2][1]-range[2][0]));
	//~ //==================2d map======================//
	//~ for(i=0;i<Cat.Len[grpid];i++)
	//~ {
		//~ pid=PIndex[i];
		//~ for(j=0;j<3;j++)
		//~ {
			//~ grid[j]=floor((Pdat.Pos[pid][j]-range[j][0])/step[j]);
			//~ if(grid[j]>=NGRID) grid[j]=NGRID-1;
		//~ }
		//~ mapxy[grid[0]][grid[1]]++;
		//~ mapxz[grid[0]][grid[2]]++;
		//~ mapyz[grid[1]][grid[2]]++;
	//~ }
	//~ 
	//~ sprintf(buf,"%s/subfindfofmapxy_%03d_%03d",outputdir,Nsnap,grpid);
	//~ if((fp=fopen(buf,"w"))==NULL)
	//~ {
		//~ fprintf(logfile,"error: file open failed for %s!\n",buf);
		//~ exit(1);
	//~ }
	//~ for(i=0;i<NGRID;i++)
	//~ {
		//~ for(j=0;j<NGRID;j++)
			//~ fprintf(fp,"%d\t",mapxy[i][j]);
		//~ fprintf(fp,"\n");
	//~ }
	//~ fclose(fp);
	
	//~ sprintf(buf,"%s/subfindfofmapxz_%03d_%03d",outputdir,Nsnap,grpid);
	//~ if((fp=fopen(buf,"w"))==NULL)
	//~ {
		//~ fprintf(logfile,"error: file open failed for %s!\n",buf);
		//~ exit(1);
	//~ }
	//~ for(i=0;i<NGRID;i++)
	//~ {
		//~ for(j=0;j<NGRID;j++)
			//~ fprintf(fp,"%d\t",mapxz[i][j]);
		//~ fprintf(fp,"\n");
	//~ }
	//~ fclose(fp);
	
	//~ sprintf(buf,"%s/subfindfofmapyz_%03d_%03d",outputdir,Nsnap,grpid);
	//~ if((fp=fopen(buf,"w"))==NULL)
	//~ {
		//~ fprintf(logfile,"error: file open failed for %s!\n",buf);
		//~ exit(1);
	//~ }
	//~ for(i=0;i<NGRID;i++)
	//~ {
		//~ for(j=0;j<NGRID;j++)
			//~ fprintf(fp,"%d\t",mapyz[i][j]);
		//~ fprintf(fp,"\n");
	//~ }
	//~ fclose(fp);
	
	//===========most bound===============//
	sprintf(buf,"%s/subfindsubcenmap_%03d_%03d",outputdir,Nsnap,grpid);
	if((fp=fopen(buf,"w"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		exit(1);
	}
	for(i=0;i<SubCat.GrpLen_Sub[grpid];i++)
	{
		pid=PIndex_all[SubMostBoundID[i]];
		for(j=0;j<3;j++)
		{
			grid[j]=floor((Pdat.Pos[pid][j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
			fprintf(fp,"%d\t",grid[j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);	
	sprintf(buf,"%s/subfindsubcen_%03d_%03d",outputdir,Nsnap,grpid);
	if((fpr=fopen(buf,"w"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		exit(1);
	}
		for(i=0;i<SubCat.GrpLen_Sub[grpid];i++)
	{
		pid=PIndex_all[SubMostBoundID[i]];
		fprintf(fpr,"%f\t%f\t%f\n",Pdat.Pos[pid][0],Pdat.Pos[pid][1],Pdat.Pos[pid][2]);
	}
	fclose(fpr);
	
	//=================min of pot==================//
	//~ sprintf(buf,"%s/subminmap_%03d_%03d",outputdir,Nsnap,grpid);
	//~ if((fp=fopen(buf,"w"))==NULL)
	//~ {
		//~ fprintf(logfile,"error: file open failed for %s!\n",buf);
		//~ exit(1);
	//~ }
	//~ sprintf(buf,"%s/submin_%03d_%03d",outputdir,Nsnap,grpid);
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
	sprintf(buf,"%s/subfindsubcommap_%03d_%03d",outputdir,Nsnap,grpid);
	if((fp=fopen(buf,"w"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		exit(1);
	}
	sprintf(buf,"%s/subfindsubcom_%03d_%03d",outputdir,Nsnap,grpid);
	if((fpr=fopen(buf,"w"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		exit(1);
	}
	for(subid=0;subid<SubCat.GrpLen_Sub[grpid];subid++)//NOTE: if not first fof, some modification for subid
	{
		////~ com[0]=0.;com[1]=0.;com[2]=0.;
		////~ for(i=0;i<SubCat.SubLen[subid];i++)
		////~ {
		//	//~ pid=SubCat.PSubArr[subid][i];
		//	//~ for(j=0;j<3;j++)
		//		//~ com[j]+=Pdat.Pos[pid][j];
		////~ }
		////~ for(j=0;j<3;j++)
		////~ {	
		//	//~ com[j]/=SubCat.SubLen[subid];
		//	//~ grid[j]=floor((com[j]-range[j][0])/step[j]);
		//	//~ if(grid[j]>=NGRID) grid[j]=NGRID-1;
		//	//~ fprintf(fp,"%d\t",grid[j]);
		// //~ }
		for(j=0;j<3;j++)
		{	
			com[j]=SubCat.CoM[subid][j];
			grid[j]=floor((com[j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
			fprintf(fp,"%d\t",grid[j]);
		}
		fprintf(fp,"\n");
		fprintf(fpr,"%d\t%f\t%f\t%f\n",SubCat.SubLen[subid],com[0],com[1],com[2]);
	}
	fclose(fp);
	fclose(fpr);
				
return 0;
}
		
