//to snapshot a subhalo with decomposition into inner part and outer part (reaccreted particles)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define NGRID 1500
	int mapxy[NGRID][NGRID]={0},mapxz[NGRID][NGRID]={0},mapyz[NGRID][NGRID]={0};

int main(int argc,char **argv)
{
	SUBCATALOGUE SubCatMin,SubCatMax;
	int SnapMin=58,SnapMax=63,i,j,pid;
	int SIDMin=7,SIDMax=6;
	FILE *fp;
	char buf[1024];
	float step[3],range[3][2];

	int *PIndexOut,NOut,*PIndexIn,NIn,*PIndexAll,NAll,grid[3];
	short *Mask;
	double *pot,*Kin,s[3],v[3];
	double Hz,sqa,Time,PartMass;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/image",SUBCAT_DIR);
	mkdir(outputdir,0755);
	logfile=stdout;

//prepare particles
	load_sub_catalogue(SnapMax,&SubCatMax,SUBCAT_DIR);
	load_sub_catalogue(SnapMin,&SubCatMin,SUBCAT_DIR);
	load_particle_data(SnapMin,SNAPSHOT_DIR);
	fresh_ID2Index(&SubCatMax,-2);	fresh_ID2Index(&SubCatMin,-2);
	Mask=calloc(NP_DM,sizeof(short));
	for(i=0;i<SubCatMin.SubLen[SIDMin];i++)
	{
		pid=SubCatMin.PSubArr[SIDMin][i];
		Mask[pid]=1;
	}
	PIndexOut=mymalloc(sizeof(int)*SubCatMax.SubLen[SIDMax]);
	NOut=0;
	for(i=0;i<SubCatMax.SubLen[SIDMax];i++)
	{
		pid=SubCatMax.PSubArr[SIDMax][i];
		if(Mask[pid]==0)
		{
			PIndexOut[NOut]=pid;
			NOut++;
		}
	}
	free(Mask);
	PIndexIn=SubCatMin.PSubArr[SIDMin];
	NIn=SubCatMin.SubLen[SIDMin];
	NAll=NOut+NIn;
	PIndexAll=mymalloc(sizeof(int)*NAll);
	memcpy(PIndexAll,PIndexIn,sizeof(int)*NIn);
	memcpy(PIndexAll+NIn,PIndexOut,sizeof(int)*NOut);

//enclosing box	
	for(i=0;i<3;i++)
		for(j=0;j<2;j++)
			range[i][j]=Pdat.Pos[PIndexOut[0]][i];
	for(i=1;i<NOut;i++)
	{
		pid=PIndexOut[i];
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
	sprintf(buf,"%s/outsize_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=0;i<3;i++)	
		fprintf(fp,"%f,%f\n",range[i][0],range[i][1]);
	fclose(fp);
	
	//~ //==================2d map======================//
	for(i=0;i<NOut;i++)
	{
		pid=PIndexOut[i];
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
	
	sprintf(buf,"%s/outmapxy_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxy[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	sprintf(buf,"%s/outmapxz_%03d_%d",outputdir,SnapMax,SIDMax);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	sprintf(buf,"%s/outmapyz_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapyz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
//==bound map==//
	for(i=0;i<NGRID;i++)
		for(j=0;j<NGRID;j++)
		{
			mapxy[i][j]=0;
			mapyz[i][j]=0;
			mapxz[i][j]=0;
		}
	for(i=0;i<NIn;i++)
	{
		pid=PIndexIn[i];
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
	
	sprintf(buf,"%s/inmapxy_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxy[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	sprintf(buf,"%s/inmapxz_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	sprintf(buf,"%s/inmapyz_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapyz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	//=================position and energy of particles==================//

	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	 sqa = sqrt(Time);
	 #endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	 for(i=0;i<3;i++)
	 {
		 s[i]=SubCatMin.Property[SIDMin].CoM[i];
		 v[i]=SubCatMin.Property[SIDMin].VCoM[i];
	 }
	 pot=mymalloc(sizeof(double)*NAll);
	 Kin=mymalloc(sizeof(double)*NAll);
	 
	//Kinetic energy
	double dx,dy,dz,dvx,dvy,dvz;
	for(i=0;i<NAll;i++)
	{
		 dvx=sqa*(Pdat.Vel[PIndexAll[i]][0])-v[0];//relative vel.
		 dvy=sqa*(Pdat.Vel[PIndexAll[i]][1])-v[1];
		 dvz=sqa*(Pdat.Vel[PIndexAll[i]][2])-v[2];
		 dx=Pdat.Pos[PIndexAll[i]][0]-s[0];
		 dy=Pdat.Pos[PIndexAll[i]][1]-s[1];
		 dz=Pdat.Pos[PIndexAll[i]][2]-s[2];
		 #ifdef PERIODIC_BDR
		 dx=NEAREST(dx);
		 dy=NEAREST(dy);
		 dz=NEAREST(dz);
		 #endif
		 dx*=Time;dy*=Time;dz*=Time;
		 dvx+=Hz*dx;//add Hubble flow
		 dvy+=Hz*dy;
		 dvz+=Hz*dz;
		 Kin[i]=0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
	 }
		 
	//potential energy
		tree_tree_allocate(TREE_ALLOC_FACTOR*NAll,NAll);
		//~ maketree(NAll,PIndexAll,Pdat.Pos);
		//~ for(i=0;i<NAll;i++)
		//~ {
		//~ pot[i]=tree_treeevaluate_potential(Pdat.Pos[PIndexAll[i]],PIndexAll,Pdat.Pos);
		//~ pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;/*exclude self-potential 
																	//~ *which was included when evaluating potential
																	//~ *  (-2.8M/h=-M/softening when r=0)*/
		//~ }
		
		maketree(NIn,PIndexIn,Pdat.Pos);
		for(i=0;i<NIn;i++)
		{
		pot[i]=tree_treeevaluate_potential(Pdat.Pos[PIndexAll[i]],PIndexIn,Pdat.Pos);
		pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;/*exclude self-potential 
																	*which was included when evaluating potential
																	*  (-2.8M/h=-M/softening when r=0)*/
		}
		for(i=NIn;i<NAll;i++)
		{
		pot[i]=tree_treeevaluate_potential(Pdat.Pos[PIndexAll[i]],PIndexIn,Pdat.Pos);
		pot[i]=(PartMass*pot[i])*G/Time;/*exclude self-potential 
										*which was included when evaluating potential
											*  (-2.8M/h=-M/softening when r=0)*/
		}
		
		tree_tree_free();
		
		
	sprintf(buf,"%s/inpart_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=0;i<NIn;i++)
	{
		for(j=0;j<3;j++)
		fprintf(fp,"%g,",Pdat.Pos[PIndexAll[i]][j]);
		fprintf(fp,"%g,%g\n",Kin[i],pot[i]);
	}	
	fclose(fp);
	sprintf(buf,"%s/outpart_%03d_%d",outputdir,SnapMin,SIDMin);
	myfopen(fp,buf,"w");
	for(i=NIn;i<NAll;i++)
	{
		for(j=0;j<3;j++)
		fprintf(fp,"%g,",Pdat.Pos[PIndexAll[i]][j]);
		fprintf(fp,"%g,%g\n",Kin[i],pot[i]);
	}	
	fclose(fp);
	
	free(Kin);
	free(pot);
				
return 0;
}
