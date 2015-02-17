//to make snapshots of subhalo particles in phase space: (v,r) and (K,P)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define SUBFIND_DIR "/home/kambrain/data/6702DM/subcatS"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat;
	int Nsnap,i,j,pid;
	int subid;
	FILE *fp;
	char buf[1024];

	int *PIndex,N;
	double *pot,*Kin,s[3],vs[3];
	double Hz,sqa,Time,PartMass;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/image",SUBCAT_DIR);
	mkdir(outputdir,0755);
	logfile=stdout;
	
	Nsnap=atoi(argv[1]);
	subid=atoi(argv[2]);
	
//prepare particles
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	
	PIndex=SubCat.PSubArr[subid];
	N=SubCat.SubLen[subid];
	pot=mymalloc(sizeof(double)*N);
	Kin=mymalloc(sizeof(double)*N);
	
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
			 vs[i]=SubCat.Property[subid].VCoM[i];
		 }
		 
		//Kinetic energy
		double dx,dy,dz,dvx,dvy,dvz;
		for(i=0;i<N;i++)
		{
			 dvx=sqa*(Pdat.Vel[PIndex[i]][0])-vs[0];//relative vel.
			 dvy=sqa*(Pdat.Vel[PIndex[i]][1])-vs[1];
			 dvz=sqa*(Pdat.Vel[PIndex[i]][2])-vs[2];
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
				
		sprintf(buf,"%s/energymap_%03d_%d.%d",outputdir,Nsnap,subid,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
			fprintf(fp,"%g,%g\n",Kin[i],pot[i]);
		fclose(fp);
		sprintf(buf,"%s/relposmap_%03d_%d.%d",outputdir,Nsnap,subid,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				fprintf(fp,"%g,",Pdat.Pos[PIndex[i]][j]-s[j]);
			fprintf(fp,"\n");
		}
		fclose(fp);
		sprintf(buf,"%s/relvelmap_%03d_%d.%d",outputdir,Nsnap,subid,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				fprintf(fp,"%g,",sqa*(Pdat.Vel[PIndex[i]][j])-vs[j]);
			fprintf(fp,"\n");
		}
		fclose(fp);
		sprintf(buf,"%s/posmap_%03d_%d.%d",outputdir,Nsnap,subid,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				fprintf(fp,"%g,",Pdat.Pos[PIndex[i]][j]);
			fprintf(fp,"\n");
		}
		fclose(fp);
		sprintf(buf,"%s/velmap_%03d_%d.%d",outputdir,Nsnap,subid,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				fprintf(fp,"%g,",sqa*(Pdat.Vel[PIndex[i]][j]));
			fprintf(fp,"\n");
		}
		fclose(fp);
		printf("CoM:\n%g,%g,%g\n",s[0],s[1],s[2]);
		printf("VCoM:\n%g,%g,%g\n",vs[0],vs[1],vs[2]);
		
		free_sub_catalogue(&SubCat);
		free(Kin);
		free(pot);
				
return 0;
}
