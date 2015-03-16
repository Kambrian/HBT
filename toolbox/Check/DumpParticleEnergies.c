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
extern void calc_particle_energy(HBTInt subid, HBTReal *Energy);

int main(int argc, char** argv)
{
	HBTInt Nsnap;
	HBTInt subid,i,j,pid;
	
	logfile=stdout;//redirect BT routines' log info to standard output
	
	Nsnap=MaxSnap-1;
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap="HBTIFMT"\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&SubCat,-2);
	free_PIDHash();
	
	printf("Nsnap=%d, Nsubs=%d\n",Nsnap, (int)SubCat.Nsubs);
	
	FILE *fp;
	myfopen(fp,"BindingEnergies.txt","w");

#ifdef HALO_PARA
omp_set_nested(0);
#endif

	HBTReal *Energy;
	Energy=mymalloc(sizeof(HBTReal)*SubCat.Nids);
#ifdef HALO_PARA
	#pragma omp parallel for
#endif
	for(subid=0;subid<SubCat.Nsubs;subid++)
	  calc_particle_energy(subid,Energy+SubCat.SubOffset[subid]);
	
	fprintf(fp, "#PID E x y z vx vy vz\n");
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
	  for(i=0;i<SubCat.SubLen[subid];i++)
	  {
	    pid=SubCat.PSubArr[subid][i];
	    fprintf(fp, HBTIFMT" %e %f %f %f %f %f %f\n",
             Pdat.PID[pid],
             Energy[SubCat.SubOffset[subid]+i],       // [10^10Msun/h (km/sec)^2], depending on unit in parameter file.
             Pdat.Pos[pid][0],                    // [kpc/h], depending on parameter file
             Pdat.Pos[pid][1],                    // [kpc/h]
             Pdat.Pos[pid][2],                    // [kpc/h]
             Pdat.Vel[pid][0],                   // [km/sec], depending on parameter file
             Pdat.Vel[pid][1],                   // [km/sec]
             Pdat.Vel[pid][2]                    // [km/sec]
             );
	  }
	}
	fclose(fp);
	myfree(Energy);
	erase_sub_catalogue(&SubCat);
	free_particle_data();
	
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
  
  tree_tree_allocate(TREE_ALLOC_FACTOR*SubCat.SubLen[subid],SubCat.SubLen[subid]);
  maketree(SubCat.SubLen[subid],SubCat.PSubArr[subid],Pdat.Pos);
  
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