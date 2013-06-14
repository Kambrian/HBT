//to track the energy evolution of subhalos from SnapMin
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat,SubCatMin;
	int Nsnap,SnapMin=43,SnapMax=99,i,j,pid;
	int SIDMin=5562,subid;
	int Nsubs,*pro2dest;
	FILE *fp;
	char buf[1024];

	int *PIndex,N;
	double *pot,*Kin,s[3],v[3];
	double Hz,sqa,Time,PartMass;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/image",SUBCAT_DIR);
	mkdir(outputdir,0755);
	logfile=stdout;

//prepare particles
	load_sub_catalogue(SnapMin,&SubCatMin,SUBCAT_DIR);
	
	PIndex=SubCatMin.PSubArr[SIDMin];
	N=SubCatMin.SubLen[SIDMin];
	subid=SIDMin;
	pot=mymalloc(sizeof(double)*N);
	Kin=mymalloc(sizeof(double)*N);
	for(Nsnap=SnapMin;Nsnap<=SnapMax;Nsnap++)
	{
		load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
		load_particle_data(Nsnap,SNAPSHOT_DIR);
		fresh_ID2Index(PIndex,N);

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
			 s[i]=SubCat.Property[subid].CoM[i];
			 v[i]=SubCat.Property[subid].VCoM[i];
		 }
		 
		//Kinetic energy
		double dx,dy,dz,dvx,dvy,dvz;
		for(i=0;i<N;i++)
		{
			 dvx=sqa*(Pdat.Vel[PIndex[i]][0])-v[0];//relative vel.
			 dvy=sqa*(Pdat.Vel[PIndex[i]][1])-v[1];
			 dvz=sqa*(Pdat.Vel[PIndex[i]][2])-v[2];
			 dx=Pdat.Pos[PIndex[i]][0]-s[0];
			 dy=Pdat.Pos[PIndex[i]][1]-s[1];
			 dz=Pdat.Pos[PIndex[i]][2]-s[2];
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
			tree_tree_allocate(TREE_ALLOC_FACTOR*N,N);
			maketree(N,PIndex,Pdat.Pos);
			for(i=0;i<N;i++)
			{
			pot[i]=tree_treeevaluate_potential(Pdat.Pos[PIndex[i]],PIndex,Pdat.Pos);
			pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;/*exclude self-potential 
																		*which was included when evaluating potential
																		*  (-2.8M/h=-M/softening when r=0)*/
			}
			tree_tree_free();
				
		sprintf(buf,"%s/energymap_%03d_%d.%d",outputdir,SnapMin,SIDMin,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
			fprintf(fp,"%g,%g\n",Kin[i],pot[i]);
		fclose(fp);
		sprintf(buf,"%s/posmap_%03d_%d.%d",outputdir,SnapMin,SIDMin,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				fprintf(fp,"%g,",Pdat.Pos[PIndex[i]][j]);
			fprintf(fp,"\n");
		}
		fclose(fp);
		
		for(i=0;i<N;i++)
		#ifdef PID_ORDERED
		PIndex[i]++;
		#else
		PIndex[i]=Pdat.PID[PIndex[i]];
		#endif
		if(Nsnap<MaxSnap-1)
		{
		load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);
		subid=pro2dest[subid];
		}
		free_pro2dest(pro2dest);
		free_sub_table(&SubCat);
	}
		free(Kin);
		free(pot);
				
return 0;
}
