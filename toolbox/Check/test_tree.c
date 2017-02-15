/*to calculate the binding energy for each particle in a sub*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

SUBCATALOGUE SubCat;

int main(int argc, char** argv)
{
	HBTInt Nsnap=32;
	HBTInt subid=22;
#define NumPart 401
	HBTInt PIds[NumPart];
	
	logfile=stdout;//redirect BT routines' log info to standard output

	char filename[1024];
	sprintf(filename, "/gpfs/data/jvbq85/HBT/data/AqA5/subcat2/postproc/Subhalo_32.22"); 
	FILE *fp;
	myfopen(fp, filename, "r");
	fread(PIds, sizeof(HBTInt), NumPart, fp);
	fclose(fp);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(PIds,NumPart);
	free_PIDHash();
	
	sprintf(filename, "/gpfs/data/jvbq85/HBT/data/AqA5/subcat2/postproc/Subhalo_32.22.pos"); 
	myfopen(fp, filename, "w");
	int i;
	for(i=0;i<NumPart;i++)
	  fprintf(fp, "%d, %g,%g,%g\n", Pdat.PID[PIds[i]], Pdat.Pos[PIds[i]][0], Pdat.Pos[PIds[i]][1],Pdat.Pos[PIds[i]][2]);
	fclose(fp);
	
	tree_tree_allocate(1*NumPart,NumPart);
	maketree(NumPart,PIds,Pdat.Pos);
	double pot=tree_treeevaluate_potential(Pdat.Pos[PIds[0]],PIds,Pdat.Pos);
	
	tree_tree_free();
	return 0;
}

void calc_particle_energy(HBTInt subid, HBTReal *Energy)
{
  if(!SubCat.SubLen[subid]) return;
  
  double Time, sqa, Hz, PartMass, s[3],v[3];
  HBTInt i;
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
  
  
  
#pragma omp parallel for
  for(i=0;i<SubCat.SubLen[subid];i++)
  {
     double dx,dy,dz,dvx,dvy,dvz,pot;
     HBTInt pid;
    pid=SubCat.PSubArr[subid][i];
    dvx=sqa*(Pdat.Vel[pid][0])-v[0];//relative vel.
    dvy=sqa*(Pdat.Vel[pid][1])-v[1];
    dvz=sqa*(Pdat.Vel[pid][2])-v[2];
    dx=Pdat.Pos[pid][0]-s[0];
    dy=Pdat.Pos[pid][1]-s[1];
    dz=Pdat.Pos[pid][2]-s[2];
    #ifdef PERIODIC_BDR
    dx=NEAREST(dx);
    dy=NEAREST(dy);
    dz=NEAREST(dz);
    #endif
    dx*=Time;dy*=Time;dz*=Time;
    dvx+=Hz*dx;//add Hubble flow
    dvy+=Hz*dy;
    dvz+=Hz*dz;
    Energy[i]=0.5*(dvx*dvx+dvy*dvy+dvz*dvz);//Kinetic energy
    pot=tree_treeevaluate_potential(Pdat.Pos[pid],SubCat.PSubArr[subid],Pdat.Pos);
    pot=(PartMass*pot+PartMass/SofteningHalo)*G/Time;//*exclude self-potential
    Energy[i]+=pot;
  }
  
  tree_tree_free();
}