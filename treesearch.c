#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

void tree_count_bin(HBTReal cen[3],HBTReal *edges,HBTInt nbin,HBTInt* bin_count,HBTInt *PIndex)
{	/*fill bin_count if the current part or node can be used,
	*bin_count should be initialized before passing here
	*return nextnode */
  union NODE *nop = 0;
  HBTReal r[3],r0,*nodecen,nodelen_half;
  HBTInt i,bin,bin_cut,bin_enc,no;

  no=NumPart;
  while(no>=0)
  {
     if(no < NumPart)		/* single particle */
	{
		r0=distance(cen,Pdat.Pos[PIndex[no]]);
		bin=binsert_asc(edges,nbin,r0);
		if(bin>=0&&bin<nbin)
			bin_count[bin]++;
	no=Nextnode[no];
	}
     else  	/*node*/
	{
		nodecen=Nodes[no].way.s;//note this is not really the center of cube but the center of Mass!!
		nodelen_half=Nodes[no].way.len/2;
		distance_point2cube(cen,nodecen,nodelen_half,r);//so this dist is not really accurate since we used CoM rather than center of cube!
		#ifdef BIN_LIN_0  //linear bin starting from 0
		bin_cut=floor(r[0]/edges[nbin]*nbin);
		bin=floor(r[1]/edges[nbin]*nbin);
		bin_enc=floor(r[2]/edges[nbin]*nbin);
		#else
		bin_cut=binsert_asc(edges,nbin,r[0]);
		bin=binsert_asc(edges,nbin,r[1]);
		bin_enc=binsert_asc(edges,nbin,r[2]);
		#endif
		if(bin==bin_cut&&bin==bin_enc)//no need to open
		{
			if(bin>=0&&bin<nbin)
			bin_count[bin]+=Nodes[no].way.mass;
			no=Nodes[no].way.sibling;
		}
		else//open node
			no=Nodes[no].way.nextnode;
	}
  }
}

HBTReal cutting_sphere_radius_cube(HBTReal cen[3],HBTReal cubecen[3],HBTReal cubelen_half)
{	//return the radius of the outer-cutting sphere centered at cen[3] for the cube with cubecen[] and cubelen=2*len_half
	HBTInt i;
	HBTReal dx[3];
	for(i=0;i<3;i++)
	{
	dx[i]=cen[i]-cubecen[i];
	#ifdef PERIODIC_BDR
	dx[i]=NEAREST(dx[i]);
	#endif
	dx[i]=fabs(dx[i]);
	if(dx[i]<=cubelen_half)
	dx[i]=0;
	else
	dx[i]-=cubelen_half;
	}
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}
HBTReal enclosing_sphere_radius_cube(HBTReal cen[3],HBTReal cubecen[3],HBTReal cubelen_half)
{	//return the radius of the outer-cutting sphere centered at cen[3] for the cube with cubecen[] and cubelen=2*len_half
	HBTInt i;
	HBTReal dx[3];
	for(i=0;i<3;i++)
	{
	dx[i]=cen[i]-cubecen[i];
	#ifdef PERIODIC_BDR
	dx[i]=NEAREST(dx[i]);
	#endif
	dx[i]=fabs(dx[i])+cubelen_half;
	}
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}

void distance_point2cube(HBTReal cen[3],HBTReal cubecen[3],HBTReal cubelen_half,HBTReal dist[3])
{	//return the radius of the outer-cutting sphere centered at cen[3] for the cube with cubecen[] and cubelen=2*len_half
	HBTInt i;
	HBTReal dx_abs[3],dx_cut[3],dx_enc[3];
	for(i=0;i<3;i++)
	{
	dx_abs[i]=cen[i]-cubecen[i];
	#ifdef PERIODIC_BDR
	dx_abs[i]=NEAREST(dx_abs[i]);
	#endif
	dx_abs[i]=fabs(dx_abs[i]);
	if(dx_abs[i]<=cubelen_half)
	dx_cut[i]=0;
	else
	dx_cut[i]=dx_abs[i]-cubelen_half;
	dx_enc[i]=dx_abs[i]+cubelen_half;
	}
	dist[0]=sqrt(dx_cut[0]*dx_cut[0]+dx_cut[1]*dx_cut[1]+dx_cut[2]*dx_cut[2]);//cutting sphere radius,nearest dist	
	dist[1]=sqrt(dx_abs[0]*dx_abs[0]+dx_abs[1]*dx_abs[1]+dx_abs[2]*dx_abs[2]);//central dist
	dist[2]=sqrt(dx_enc[0]*dx_enc[0]+dx_enc[1]*dx_enc[1]+dx_enc[2]*dx_enc[2]);//enclosing sphere radius,farthest dis
}

HBTInt binsert_asc(HBTReal *edges,HBTInt nbin,HBTReal target)
{//binary insert into an ascending arr,return id where edges[id]=<target<edges[id+1],id=nbin when between [edges[nbin],inf);
	HBTInt ids,ide,idm;
	if(target==edges[0])
	return 0;
	if(target>=edges[nbin])
	return nbin;//after bins
	if (target<edges[0])
	return -1;//before bins
	
	ids=0;ide=nbin;
	while(ide>ids+1)//have space for division
	{
	idm=(ids+ide)/2;
	if(target==edges[idm])
	return idm;
	if(target>edges[idm])
		ids=idm;
	else
		ide=idm;
	}
	return ids;
}

/* the following neighbour-search code also modified from SUBFIND*/
static HBTReal *NgbR2;
static HBTInt *NgbID;
static HBTInt NgbNMax;
static HBTInt NgbNMax0;
#pragma omp threadprivate(NgbR2,NgbID,NgbNMax)

HBTReal guess_ngb_range(HBTInt NumNgb)
{
HBTReal hguess,comoving_num_dens;
HBTInt NumNgbSafe;
//~ if(NumNgb<5)
//~ NumNgbSafe=ceil(pow((1+sqrt(1.0+4*NumNgb))/2.0,2));//so that NumNgbSafe-sqrt(NumNgbSafe)=NumNgb,to allow for Poisson fluctuation.
//~ else
NumNgbSafe=NumNgb;

NgbNMax0=256*NumNgbSafe;
comoving_num_dens = header.Omega0 * 3 * HUBBLE0 * HUBBLE0 / (8 * 3.141593 * G) / header.mass[1];//comoving number density for DM,also for gas if NP_gas=NP_DM
printf("Dens=%g\n",comoving_num_dens);
hguess = pow(3 * NumNgbSafe / (4 * 3.141593) / comoving_num_dens, 1.0 / 3);
return hguess;
}

HBTReal guess_ngb_range_halo(HBTInt NumNgb)
{//for guess ngb_range inside halos; use 200*mean_density as background
HBTReal hguess,comoving_num_dens;
HBTInt NumNgbSafe;
//~ if(NumNgb<5)
//~ NumNgbSafe=ceil(pow((1+sqrt(1.0+4*NumNgb))/2.0,2));//so that NumNgbSafe-sqrt(NumNgbSafe)=NumNgb,to allow for Poisson fluctuation.
//~ else
NumNgbSafe=NumNgb;

NgbNMax0=64*NumNgbSafe;
comoving_num_dens = header.Omega0 * 3 * HUBBLE0 * HUBBLE0 / (8 * 3.141593 * G) / header.mass[1];//comoving number density for DM,also for gas if NP_gas=NP_DM
printf("Dens=%g\n",comoving_num_dens);
hguess = pow(3 * NumNgbSafe / (4 * 3.141593) / comoving_num_dens /200., 1.0 / 3);
return hguess;
}

HBTInt treesearch_nearest(HBTReal cen[3],HBTReal hguess,HBTInt *PIndex,HBTReal PPos[][3])
{//find the nearest dark matter neighbor to the target position cen[3]
	HBTInt numngb,pid;
	NgbNMax=NgbNMax0;
	NgbR2=mymalloc(sizeof(HBTReal)*NgbNMax);
	NgbID=mymalloc(sizeof(HBTInt)*NgbNMax);
	
	numngb = treesearch_sphere(cen, hguess, PIndex,PPos);
	while(numngb<1)
    {
	 hguess *= 1.26;//double the guess volume
     numngb = treesearch_sphere(cen, hguess, PIndex,PPos);
    }
	pid=NgbID[min_of_vec(NgbR2,numngb)];	
	free(NgbR2);
	free(NgbID);
	return pid;
}

#ifndef PERIODIC_BDR
#define NEAREST
#endif

HBTInt treesearch_sphere(HBTReal searchcenter[3], HBTReal radius,HBTInt *PIndex,HBTReal PPos[][3])
{/*find a list of particles from PIndex within radius around searchcenter
  * return the number of neighbors
  * also store distance^2 and ID in the global variables NgbR2 and NgbID */
  HBTInt numngb, no, p;
  double dx, dy, dz, r2, h2;
  union NODE *this;

  h2 = radius * radius;

  numngb = 0;
  no = NumPart;

  while(no >= 0)
    {
      if(no < NumPart)		/* single particle */
	{
	  p = PIndex[no];
	  no = Nextnode[no];

	  dx = PPos[p][0] - searchcenter[0];
	  #ifdef PERIODIC_BDR
	  dx=NEAREST(dx);
	  #endif
	  if(dx < -radius)
	    continue;
	  if(dx > radius)
	    continue;

	  dy =PPos[p][1] - searchcenter[1];
	  #ifdef PERIODIC_BDR
	  dy=NEAREST(dy);
	  #endif
	  if(dy < -radius)
	    continue;
	  if(dy > radius)
	    continue;

	  dz = PPos[p][2] - searchcenter[2];
	  #ifdef PERIODIC_BDR
	  dz=NEAREST(dz);
	  #endif
	  if(dz < -radius)
	    continue;
	  if(dz > radius)
	    continue;

	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      NgbR2[numngb]= r2;
	      NgbID[numngb]= p;
	      numngb++;
		  if(numngb==NgbNMax)
		  { 
			NgbNMax*=2;
		  	NgbR2=realloc(NgbR2,sizeof(HBTReal)*NgbNMax);
			NgbID=realloc(NgbID,sizeof(HBTInt)*NgbNMax);
		  }
	    }
	}
      else
	{
	  this = &Nodes[no];

	  no = Nodes[no].way.sibling;	/* in case the node can be discarded */
	//since way.s[3] is CoM rather than center of cube,compare with Len rather than Lenhalf to allow for misaligned CoM
	  if((NEAREST(this->way.s[0] - searchcenter[0]) + this->way.len) < -radius)
	    continue;
	  if((NEAREST(this->way.s[0] - searchcenter[0]) - this->way.len) > radius)
	    continue;
	  if((NEAREST(this->way.s[1] - searchcenter[1]) + this->way.len) < -radius)
	    continue;
	  if((NEAREST(this->way.s[1] - searchcenter[1]) - this->way.len) > radius)
	    continue;
	  if((NEAREST(this->way.s[2] - searchcenter[2]) + this->way.len) < -radius)
	    continue;
	  if((NEAREST(this->way.s[2] - searchcenter[2]) - this->way.len) > radius)
	    continue;

	  no = this->way.nextnode;	/* ok, we need to open the node */
	}
    }

  /*
     fprintf(Logfile,"numngb=%d\n", numngb);
   */
  return numngb;
}

#define SPH_DENS_NGB 64

double sph_density(HBTReal cen[3],HBTReal *p2hguess,HBTInt *PIndex,HBTReal PPos[][3])
{
  HBTReal hguess;
  HBTInt i, n;
  double h, hinv3, wk, u, r, rho;
  //~ clock_t T[10];
  
	HBTInt numngb;
	NgbNMax=NgbNMax0;
	NgbR2=mymalloc(sizeof(HBTReal)*NgbNMax);
	NgbID=mymalloc(sizeof(HBTInt)*NgbNMax);
	hguess=*p2hguess;
	
	//~ T[0]=clock();
	numngb = treesearch_sphere(cen, hguess, PIndex,PPos);
	//~ T[1]=clock();
	//~ printf("First search: %ld; %d ngbs found\n",T[1]-T[0], numngb);fflush(stdout);
	while(numngb<SPH_DENS_NGB)
    {
	 if(numngb)	
	 hguess *= pow((HBTReal)SPH_DENS_NGB/(HBTReal)numngb,1.0/3.0)*1.1;//update hguess adaptively, and conservatively to keep it slightly larger
	 else  //zero ngb, double hguess
	 hguess *= 2.;
	 
     numngb = treesearch_sphere(cen, hguess, PIndex,PPos);
	 //~ printf("N=%d,h=%f\n",numngb,hguess);fflush(stdout);
    }
	*p2hguess=hguess*powf((HBTReal)SPH_DENS_NGB/(HBTReal)numngb,1.0/3.0)*1.1;//to return a slight larger best guess
	//~ T[2]=clock();
	//~ printf("Search done: %ld, hguess=%f\n",T[2]-T[1],hguess);fflush(stdout);
	h=psort(SPH_DENS_NGB,numngb,NgbR2);//NgbR2 has now been partly sorted,with respect to h	
	//~ T[3]=clock();
	//~ printf("NgbSorted: %ld\n",T[3]-T[2]);fflush(stdout);
	h=sqrtf(h);
	hinv3 = 1.0 / (h * h * h);

	  for(n = 0, rho = 0; n < SPH_DENS_NGB; n++)
	    {
	      r = sqrtf(NgbR2[n]);
	      u = r / h;

	      if(u < 0.5)
		wk = hinv3 * (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
	      else
		wk = hinv3 * 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

	      rho += wk;
	    }
	free(NgbR2);
	free(NgbID);
	//~ T[4]=clock();
	//~ printf("DensCalc: %ld\n",T[4]-T[3]);fflush(stdout);
	//~ T[5]=T[4]-T[0];
	//~ printf("Summary: %f, %f, %f, %f\n",(HBTReal)(T[1]-T[0])/(HBTReal)T[5],(HBTReal)(T[2]-T[1])/(HBTReal)T[5],(HBTReal)(T[3]-T[2])/(HBTReal)T[5],(HBTReal)(T[4]-T[3])/(HBTReal)T[5]);
	//~ printf("%ld,%ld,%ld,%ld,%ld\n",T[0],T[1],T[2],T[3],T[4]);
 return rho;
}

static HBTInt *InfectionStack, *pInfected;
static double LinkLength, LinkLength2;
#pragma omp threadprivate(InfectionStack,pInfected,LinkLength,LinkLength2)

HBTInt treesearch_linkgrp(HBTReal radius, HBTInt PIndex[], struct GroupData *GrpData)
/* To link DM particles into groups
 * Input: radius: link radius
 * 		  PIndex[]: list of particles to link
 *        grpdata: pointer to GroupData structure.
 * Output: filled grpdata.
 * Return value: number of groups found, down to mass=1 (diffuse particles)
 * */
{
HBTInt i,grpid;

//==init
GrpData->GrpTags=mymalloc(sizeof(struct ParticleGroup)*GrpData->Np);
for(i=0;i<GrpData->Np;i++) 
{
	GrpData->GrpTags[i].PIndex=PIndex[i];
	GrpData->GrpTags[i].GrpID=-1;
}
GrpData->GrpLen=mymalloc(sizeof(HBTInt)*GrpData->Np);	

InfectionStack=mymalloc(sizeof(HBTInt)*GrpData->Np);//at most Np particles would be infected
pInfected=InfectionStack;//initial position
LinkLength=radius;
LinkLength2=radius*radius;

fprintf(logfile,"constructing tree...\n");fflush(logfile);
tree_tree_allocate(TREE_ALLOC_FACTOR*(size_t)GrpData->Np,GrpData->Np);
maketree(GrpData->Np,PIndex,Pdat.Pos);
	
fprintf(logfile,"Linking Groups...\n");fflush(logfile);
//~ printf("P%08d G%08d",0,0);
grpid=0;
HBTInt j,jj=1; j=NumPart/100;
fprintf(logfile,"00%%");fflush(logfile);
for(i=0;i<NumPart;i++)
{
	if(GrpData->GrpTags[i].GrpID<0)
	{
		//~ int j;
		//~ for(j=0;j<19;j++) printf("\b");
		//~ printf("P%08d G%08d",i,grpid);fflush(stdout);
		if(i>=j){fprintf(logfile,"\b\b\b%02d%%",(int)jj);fflush(logfile);jj++;j=(NumPart/100.)*jj;}
		GrpData->GrpTags[i].GrpID=grpid; //infect the seed particle
		GrpData->GrpLen[grpid]=1+treesearch_infect_particles(i,grpid, GrpData->GrpTags, Pdat.Pos);
		grpid++;
		
	}
}

tree_tree_free();

GrpData->Ngrp=grpid;
GrpData->GrpLen=realloc(GrpData->GrpLen,sizeof(HBTInt)*grpid);
fprintf(logfile,"Found "HBTIFMT" Groups\n",grpid);fflush(logfile);

if(pInfected!=InfectionStack)
{
	fprintf(logfile,"Error: InfectionStack not cleaned in treesearch_linkgrp, current pos=%d\n",
						(int)(pInfected-InfectionStack));
}
myfree(InfectionStack);

return grpid;	
}

HBTInt treesearch_infect_particles(HBTInt seed, HBTInt grpid,
		struct ParticleGroup *GrpTags, HBTReal PPos[][3])
{
/*tag all the particles that are linked to seed with grpid
 * Note if system stack size is too small, this recursive routine may crash
 * in that case you should set:  ulimit -s unlimited  (bash) before running.
**/
 HBTInt numngb, totnumngb, no, p;
  double dx, dy, dz, r2, h2;
  union NODE *this;
  HBTReal *searchcenter;

  searchcenter=PPos[GrpTags[seed].PIndex];
  //~ h2 = radius * radius;

  numngb = 0; 
  no = NumPart;

  while(no >= 0)//infect neighbours
    {
		
      if(no < NumPart)		/* single particle */
	  {
		  *pInfected = no;
		  no = Nextnode[*pInfected];
		  if(GrpTags[*pInfected].GrpID>=0) //already tagged
			continue;  
		  
		  p = GrpTags[*pInfected].PIndex;

		  dx = PPos[p][0] - searchcenter[0];
		  #ifdef PERIODIC_BDR
		  dx=NEAREST(dx);
		  #endif
		  if(dx < -LinkLength)
			continue;
		  if(dx > LinkLength)
			continue;

		  dy =PPos[p][1] - searchcenter[1];
		  #ifdef PERIODIC_BDR
		  dy=NEAREST(dy);
		  #endif
		  if(dy < -LinkLength)
			continue;
		  if(dy > LinkLength)
			continue;

		  dz = PPos[p][2] - searchcenter[2];
		  #ifdef PERIODIC_BDR
		  dz=NEAREST(dz);
		  #endif
		  if(dz < -LinkLength)
			continue;
		  if(dz > LinkLength)
			continue;

		  r2 = dx * dx + dy * dy + dz * dz;

		  if(r2 < LinkLength2) //confirm the infection (fill the stack)
			{
			  GrpTags[*pInfected].GrpID=grpid;
			  numngb++;
			  pInfected++;
			}
	  }
      else
	  {
	  this = &Nodes[no];

	  no = Nodes[no].way.sibling;	/* in case the node can be discarded */
	//since way.s[3] is CoM rather than center of cube,compare with Len rather than Lenhalf to allow for misaligned CoM
	  if((NEAREST(this->way.s[0] - searchcenter[0]) + this->way.len) < -LinkLength)
	    continue;
	  if((NEAREST(this->way.s[0] - searchcenter[0]) - this->way.len) > LinkLength)
	    continue;
	  if((NEAREST(this->way.s[1] - searchcenter[1]) + this->way.len) < -LinkLength)
	    continue;
	  if((NEAREST(this->way.s[1] - searchcenter[1]) - this->way.len) > LinkLength)
	    continue;
	  if((NEAREST(this->way.s[2] - searchcenter[2]) + this->way.len) < -LinkLength)
	    continue;
	  if((NEAREST(this->way.s[2] - searchcenter[2]) - this->way.len) > LinkLength)
	    continue;

	  no = this->way.nextnode;	/* ok, we need to open the node */
	  }
    }
	
	//~ printf("ngb=%d ",numngb);
	totnumngb=numngb;
	while(numngb>0)//pop the stack
	{
		pInfected--;
		totnumngb+=treesearch_infect_particles(*pInfected,grpid, GrpTags, PPos);
		numngb--;
	}
	return totnumngb; //total number of infected particles (excluding the seed particle)
}
