#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

static int comp_erg(const void *a, const void *b)//used to sort energy in ascending order; note that the most bound will come first( energy <0)
{
  if(((struct Energy *) a)->Erg > ((struct Energy *) b)->Erg)
    return +1;

  if(((struct Energy *) a)->Erg < ((struct Energy *) b)->Erg)
    return -1;

  return 0;
}

HBTInt unbind(HBTInt *P2Len,HBTInt **P2PIndex, struct SubProperty *Prop,HBTInt *P2Len_removed, HBTInt **P2PIndex_removed,HBTReal CoreFrac) /*P2Len=&Len, *P2PIndex=PIndex, 
																							*where PIndex[Len] is the array of size Len; 
																							* both will be updated after unbind
																							* *P2PIndex_removed=PIndex_removed need 
																							* not to be allocated as input*/
{
	HBTInt * TmpPIndex, *TmpI, *SubPIndex,*SubPIndex_removed;
	HBTReal * TmpBindingE;
	struct Energy *SubEsort;
	HBTInt i,Nbound,SubLen,SubLen_removed,CoreLen;
	double Hz,sqa,Time,PartMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz,*pot,E,Erelax;
	double SubPot,SubKin,SubAMx,SubAMy,SubAMz;
	#ifdef PERIODIC_BDR
	double mx,my,mz;//position of minimum potential
	#endif
	//time_t time_start,time_end,t1,t2,t3,t4,tt1,tt2,tt3,tt4,ttt2,ttt3; //for program timing
	//~ FILE *fp;
	//fp=fopen("converg.dat","w");
	if(*P2Len<NBOUNDMIN)
	{
		*P2PIndex_removed=*P2PIndex;//this is necessary to make legal the free() of PIndex_removed outside.
		*P2Len_removed=*P2Len;
		*P2PIndex=NULL;
		*P2Len=0;
		Prop->CoM[0]=Prop->CoM[1]=Prop->CoM[2]=0.;
		Prop->VCoM[0]=Prop->VCoM[1]=Prop->VCoM[2]=0.;
		Prop->Pot=Prop->Kin=0.;
		Prop->AM[0]=Prop->AM[1]=Prop->AM[2]=0.;
		return 0;
	}
	
	//time_start=time(NULL);
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	 sqa = sqrt(Time);
	 #endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	Nbound=*P2Len;
	SubPIndex=*P2PIndex;
	SubLen_removed=0;
	SubPIndex_removed=mymalloc(sizeof(HBTInt)*Nbound);
	
	//fprintf(fp,"%d\n",Nbound);
	//if(Nbound>NParaMin) printf("%d\n",Nbound);
	SubLen=Nbound*2+1;//just to start the loop
	TmpBindingE=mymalloc(1);//just to start the loop
	TmpI=mymalloc(1);
	pot=mymalloc(1);
	while(Nbound<SubLen*PrecMass)
	{	
		free(pot);
		free(TmpBindingE);
		free(TmpI);
		//fprintf(logfile,"Nbound=%d\n",Nbound);fflush(stdout);
		SubLen=Nbound;
		Nbound=0;
		TmpPIndex=mymalloc(sizeof(HBTInt)*SubLen);//this has been taken over by SubPIndex
		TmpI=mymalloc(sizeof(HBTInt)*SubLen);
		TmpBindingE=mymalloc(sizeof(HBTReal)*SubLen);
		pot=mymalloc(sizeof(double)*SubLen);
		
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubPIndex,Pdat.Pos);
		vx=vy=vz=sx=sy=sz=0.;
	
		#pragma omp parallel  if(SubLen>NParaMin) 
		{//\\start para
			#pragma omp for schedule(dynamic,1)
			for(i=0;i<SubLen;i++)
			{
				pot[i]=tree_treeevaluate_potential(Pdat.Pos[SubPIndex[i]],SubPIndex,Pdat.Pos);
				pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;/*exclude self-potential 
																		*which was included when evaluating potential
																		*  (-2.8M/h=-M/softening when r=0)*/
			}
			#pragma omp single
			{
			SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
			}
			#pragma omp for schedule(dynamic,1)
			for(i=0;i<SubLen;i++)
			{
			SubEsort[i].PID=SubPIndex[i];
			SubEsort[i].Erg=pot[i];
			}
			#pragma omp single
			{
			qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
			CoreLen=SubLen*CoreFrac;
			CoreLen=(CoreLen>CoreLenMin)?CoreLen:CoreLenMin;
			#ifdef PERIODIC_BDR
			mx=Pdat.Pos[SubEsort[0].PID][0];
			my=Pdat.Pos[SubEsort[0].PID][1];
			mz=Pdat.Pos[SubEsort[0].PID][2];
			#endif
			//~ printf("corelen %d\n",CoreLen);
			}
			#pragma omp for private(i,dx,dy,dz) reduction (+:vx,vy,vz,sx,sy,sz) schedule(dynamic,1)
			for(i=0;i<CoreLen;i++)
			{
			vx+=Pdat.Vel[SubEsort[i].PID][0];
			vy+=Pdat.Vel[SubEsort[i].PID][1];
			vz+=Pdat.Vel[SubEsort[i].PID][2];
			#ifdef PERIODIC_BDR
			dx=Pdat.Pos[SubEsort[i].PID][0]-mx;
			dy=Pdat.Pos[SubEsort[i].PID][1]-my;
			dz=Pdat.Pos[SubEsort[i].PID][2]-mz;
			sx+=NEAREST(dx);
			sy+=NEAREST(dy);
			sz+=NEAREST(dz);
			#else
			sx+=Pdat.Pos[SubEsort[i].PID][0];
			sy+=Pdat.Pos[SubEsort[i].PID][1];
			sz+=Pdat.Pos[SubEsort[i].PID][2];
			#endif
			}
			#pragma omp single
			{
			vx/=CoreLen;
			vy/=CoreLen;
			vz/=CoreLen;
			sx/=CoreLen;
			sy/=CoreLen;
			sz/=CoreLen;
			free(SubEsort);
			#ifdef PERIODIC_BDR
			sx+=mx;
			sy+=my;
			sz+=mz;
			#endif
			}
			#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic,1) 
			for(i=0;i<SubLen;i++)
			{
			 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-vx);//relative vel.
			 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-vy);
			 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-vz);
			 dx=Pdat.Pos[SubPIndex[i]][0]-sx;
			 dy=Pdat.Pos[SubPIndex[i]][1]-sy;
			 dz=Pdat.Pos[SubPIndex[i]][2]-sz;
			 #ifdef PERIODIC_BDR
			 dx=NEAREST(dx);
			 dy=NEAREST(dy);
			 dz=NEAREST(dz);
			 #endif
			 dx*=Time;dy*=Time;dz*=Time;
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			 E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
			 #pragma omp critical (unbinding)
			#ifdef NO_UNBINDING
				if(1)
			#else 
				#ifdef E_Relax
				if(E<-E_Relax*pot[i]) //(E_Relax+1)*pot[i]+K<0
				#else
				if(E<0)
				#endif
			#endif
				{
				TmpI[Nbound]=i;
				TmpBindingE[Nbound]=E;
				Nbound++;
				}
				else
				{
				SubPIndex_removed[SubLen_removed]=SubPIndex[i];
				SubLen_removed++;	
				}
			}
			#pragma omp for schedule(dynamic,1) 
			for(i=0;i<Nbound;i++)
				TmpPIndex[i]=SubPIndex[TmpI[i]];
		}//\\end para
		//fprintf(fp,"%d\n",Nbound);
		tree_tree_free();
		free(SubPIndex);//am i also freeing SubCatTmp->PSubArr[subhaloid]??............................Ok,fine.
		SubPIndex=TmpPIndex;//take over the bound part to begin a new loop 
		if(Nbound<NBOUNDMIN)
		{
			memcpy(SubPIndex_removed+SubLen_removed,TmpPIndex,Nbound*sizeof(HBTInt));
			*P2PIndex_removed=SubPIndex_removed;
			*P2Len_removed=SubLen_removed+Nbound;
			free(TmpPIndex);
			free(TmpBindingE);
			free(pot);
			free(TmpI);
			*P2Len=0;
			*P2PIndex=NULL;
			Prop->CoM[0]=Prop->CoM[1]=Prop->CoM[2]=0.;
			Prop->VCoM[0]=Prop->VCoM[1]=Prop->VCoM[2]=0.;
			Prop->Pot=Prop->Kin=0.;
			Prop->AM[0]=Prop->AM[1]=Prop->AM[2]=0.;
			return 0;//no bound sub found, return 0;
		}			
	}
	//t1=time(NULL);
	
	//sort the particles and return the IDs
	SubLen=Nbound;//NOTE: this is only necessary when using tolerated (approximated) binding criteria
	SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
	for(i=0;i<SubLen;i++)
	{
		SubEsort[i].PID=SubPIndex[i];
		SubEsort[i].Erg=TmpBindingE[i];
	}
	qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
	for(i=0;i<SubLen;i++)
		SubPIndex[i]=SubEsort[i].PID;//copy the sorted index back

	free(TmpBindingE);
	free(SubEsort);
	*P2Len=SubLen;
	*P2PIndex=SubPIndex;
	
	//return properties of the bound structure
	Prop->CoM[0]=sx;Prop->CoM[1]=sy;Prop->CoM[2]=sz; 
	Prop->VCoM[0]=sqa*vx;Prop->VCoM[1]=sqa*vy;Prop->VCoM[2]=sqa*vz; 
	SubPot=SubKin=SubAMx=SubAMy=SubAMz=0.;
	#pragma omp parallel for private(i,dvx,dvy,dvz,dx,dy,dz) schedule(dynamic,1) reduction(+:SubPot,SubKin,SubAMx,SubAMy,SubAMz)
	for(i=0;i<SubLen;i++)
	{
	 SubPot+=pot[TmpI[i]];
	 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-vx);//relative vel.
	 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-vy);
	 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-vz);
	 dx=Pdat.Pos[SubPIndex[i]][0]-sx;
	 dy=Pdat.Pos[SubPIndex[i]][1]-sy;
	 dz=Pdat.Pos[SubPIndex[i]][2]-sz;
	 #ifdef PERIODIC_BDR
	 dx=NEAREST(dx);
	 dy=NEAREST(dy);
	 dz=NEAREST(dz);
	 #endif
	 dx*=Time;dy*=Time;dz*=Time;
	 dvx+=Hz*dx;//add Hubble flow
	 dvy+=Hz*dy;
	 dvz+=Hz*dz;
	 SubKin+=0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
	 SubAMx+=dy*dvz-dz*dvy;
	 SubAMy+=dx*dvz-dz*dvx;
	 SubAMz+=dx*dvy-dy*dvx;
	}
	Prop->Pot=SubPot/SubLen;
	Prop->Kin=SubKin/SubLen;
	Prop->AM[0]=SubAMx/SubLen;Prop->AM[1]=SubAMy/SubLen;Prop->AM[2]=SubAMz/SubLen;
	
	free(pot);
	free(TmpI);
	//return the removed particles
	if((*P2Len_removed=SubLen_removed))
		*P2PIndex_removed=realloc(SubPIndex_removed,sizeof(HBTInt)*SubLen_removed);
		else
		{
		free(SubPIndex_removed);
		*P2PIndex_removed=NULL;
		}
	//time_end=time(NULL);
	//printf("%ld (%f),%ld (%f) sec\n",t1-time_start,(HBTReal)(t1-time_start)/(time_end-time_start),time_end-t1,(HBTReal)(time_end-t1)/(time_end-time_start));
	return 1;//found bound structure, return 1;
}

void unbind_sub_recursive(HBTInt mainsubID,HBTInt *P2Len_removed,HBTInt **P2PIndex_removed,SUBCATALOGUE *SubCat,SRCCATALOGUE *SrcCat)
{
	HBTInt son,sib,sonN_removed,*sonPIndex_removed;
	HBTInt Len,*PIndex;
	
	if((son=SubCat->sub_hierarchy[mainsubID].sub)>=0)//has sons,unbind them recursively to update its reservior
	{	
		Len=SubCat->SubLen[mainsubID];
		PIndex=SubCat->PSubArr[mainsubID];
		sib=son;
		while(sib>=0)
		{
			unbind_sub_recursive(sib,&sonN_removed,&sonPIndex_removed,SubCat,SrcCat);
			//add the removed particles to the nibs sub
			if(sonN_removed)
			{
			PIndex=realloc(PIndex,sizeof(HBTInt)*(Len+sonN_removed));
			memcpy(PIndex+Len,sonPIndex_removed,sonN_removed*sizeof(HBTInt));
			free(sonPIndex_removed);
			Len+=sonN_removed;
			}
			sib=SubCat->sub_hierarchy[sib].next;
		}
		SubCat->SubLen[mainsubID]=Len;
		SubCat->PSubArr[mainsubID]=PIndex;

	}
	//now unbind the mainsub itself
	unbind(SubCat->SubLen+mainsubID,SubCat->PSubArr+mainsubID,SubCat->Property+mainsubID,P2Len_removed,P2PIndex_removed,SrcCat->CoreFrac[mainsubID]);
	narrow_srccat(SrcCat,SubCat, mainsubID);
}
