#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"

static int comp_erg(const void *a, const void *b)//used to sort energy in ascending order; note that the most bound will come first( energy <0)
{
  if(((struct Energy *) a)->Erg > ((struct Energy *) b)->Erg)
    return +1;

  if(((struct Energy *) a)->Erg < ((struct Energy *) b)->Erg)
    return -1;

  return 0;
}
HBTInt unbindgas(HBTInt *P2Len,HBTInt **P2PIndex, struct GasProperty *Prop, HBTInt SubLen, HBTInt *SubArr,HBTReal SubCoM[3],HBTReal SubVCoM[3]) /*P2Len=&Len, *P2PIndex=PIndex, 
																*where PIndex[Len] is the array of size Len; 
																* both will be updated after unbind*/
{
	HBTInt * GasArr,*BndGasArr, *BndGasInd;
	HBTReal * GasE, E;
	struct Energy *SubEsort;
	HBTInt i,j,pid,gaslen,gasbound;
	double Hz,sqa,Time,DMMass,GasMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz,*pot;
	double SubPot,SubKin,AMx,AMy,AMz,U;

	if((0==SubLen)||(0==(*P2Len))) //do nothing for null subs
	{
		*P2Len=0;
		*P2PIndex=NULL;
		Prop->CoM[0]=Prop->CoM[1]=Prop->CoM[2]=0.;
		Prop->VCoM[0]=Prop->VCoM[1]=Prop->VCoM[2]=0.;
		Prop->Pot=Prop->Kin=0.;
		Prop->AM[0]=Prop->AM[1]=Prop->AM[2]=0.;
		Prop->U=0.;
		return 0;
	}
	
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	Hz=header.Hz;
	GasMass=header.mass[0];
	DMMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	gaslen=*P2Len;
	GasArr=*P2PIndex;
	gasbound=0;
	BndGasArr=mymalloc(sizeof(HBTInt)*gaslen);
	BndGasInd=mymalloc(sizeof(HBTInt)*gaslen);
	GasE=mymalloc(sizeof(HBTReal)*gaslen);
	pot=mymalloc(sizeof(double)*gaslen);

	tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
	maketree(SubLen,SubArr,Pdat.Pos);

	sx=sy=sz=vx=vy=vz=0.;
	SubPot=SubKin=AMx=AMy=AMz=U=0.;
	#pragma omp parallel if(SubLen>NParaMin) 
	{
	#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic,1) 
	for(i=0;i<gaslen;i++)
	{		
		dvx=sqa*Gdat.Vel[GasArr[i]][0]-SubVCoM[0];//relative vel.Note: SubVCoM is physical Vel
		dvy=sqa*Gdat.Vel[GasArr[i]][1]-SubVCoM[1];
		dvz=sqa*Gdat.Vel[GasArr[i]][2]-SubVCoM[2];
		dx=Gdat.Pos[GasArr[i]][0]-SubCoM[0];
		dy=Gdat.Pos[GasArr[i]][1]-SubCoM[1];
		dz=Gdat.Pos[GasArr[i]][2]-SubCoM[2];
		#ifdef PERIODIC_BDR
		dx=NEAREST(dx);
		dy=NEAREST(dy);
		dz=NEAREST(dz);
		#endif
		dx*=Time;dy*=Time;dz*=Time;
		dvx+=Hz*dx;//add Hubble flow
		dvy+=Hz*dy;
		dvz+=Hz*dz;
		
		pot[i]=tree_treeevaluate_potential(Gdat.Pos[GasArr[i]],SubArr,Pdat.Pos);
		pot[i]=DMMass*pot[i]*G/Time;//no self energy since this is gas in dm pot
		
		#ifdef THERMAL_BOUND	
		E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz)+Gdat.U[GasArr[i]];//necessary to add internal energy?
		#else
		E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);//necessary to add internal energy?
		#endif
		
		#pragma omp critical (unbinding) 
		if(E<0)
		{
		BndGasInd[gasbound]=i;
		GasE[gasbound]=E;
		gasbound++;
		}
	}
	#pragma omp for private(i,pid,dvx,dvy,dvz,dx,dy,dz) \
							reduction(+:sx,sy,sz,vx,vy,vz,SubPot,SubKin,AMx,AMy,AMz,U) \
							schedule(dynamic,1) 
	for(i=0;i<gasbound;i++)
	{
		SubPot+=pot[BndGasInd[i]];
		
		pid=GasArr[BndGasInd[i]];
		BndGasArr[i]=pid;
		sx+=Gdat.Pos[pid][0];
		sy+=Gdat.Pos[pid][1];
		sz+=Gdat.Pos[pid][2];
		vx+=Gdat.Vel[pid][0];
		vy+=Gdat.Vel[pid][1];
		vz+=Gdat.Vel[pid][2];
		U+=Gdat.U[pid];
				
		dvx=sqa*Gdat.Vel[pid][0]-SubVCoM[0];//relative vel.Note: VCoM is physical Vel
		dvy=sqa*Gdat.Vel[pid][1]-SubVCoM[1];
		dvz=sqa*Gdat.Vel[pid][2]-SubVCoM[2];
		dx=Gdat.Pos[pid][0]-SubCoM[0];
		dy=Gdat.Pos[pid][1]-SubCoM[1];
		dz=Gdat.Pos[pid][2]-SubCoM[2];
		#ifdef PERIODIC_BDR
		dx=NEAREST(dx);
		dy=NEAREST(dy);
		dz=NEAREST(dz);
		#endif
		dx*=Time;dy*=Time;dz*=Time;
		dvx+=Hz*dx;//add Hubble flow
		dvy+=Hz*dy;
		dvz+=Hz*dz;
		
		AMx+=dy*dvz-dz*dvy;
		AMy+=dx*dvz-dz*dvx;
		AMz+=dx*dvy-dy*dvx;
		SubKin+=0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
	}
	}
	if(gasbound)
	{
	Prop->CoM[0]=sx/gasbound;Prop->CoM[1]=sy/gasbound;Prop->CoM[2]=sz/gasbound;
	Prop->VCoM[0]=sqa*vx/gasbound;Prop->VCoM[1]=sqa*vy/gasbound;Prop->VCoM[2]=sqa*vz/gasbound;
	Prop->Pot=SubPot/gasbound;
	Prop->Kin=SubKin/gasbound;
	Prop->AM[0]=AMx/gasbound;Prop->AM[1]=AMy/gasbound;Prop->AM[2]=AMz/gasbound;
	Prop->U=U/gasbound;
	}
	else
	{
	Prop->CoM[0]=Prop->CoM[1]=Prop->CoM[2]=0.;
	Prop->VCoM[0]=Prop->VCoM[1]=Prop->VCoM[2]=0.;
	Prop->Pot=0.;
	Prop->Kin=0.;
	Prop->AM[0]=Prop->AM[1]=Prop->AM[2]=0.;
	Prop->U=0.;
	}
	tree_tree_free();
	free(GasArr);
	free(BndGasInd);
	free(pot);
	if(gasbound)
	GasArr=realloc(BndGasArr,sizeof(HBTInt)*gasbound);//take over the bound part to begin a new loop 
	else
	GasArr=NULL;
	gaslen=gasbound;
		
	SubEsort=mymalloc(sizeof(struct Energy)*gaslen);
	for(i=0;i<gaslen;i++)
	{
		SubEsort[i].PID=GasArr[i];
		SubEsort[i].Erg=GasE[i];
	}
	qsort(SubEsort,gaslen,sizeof(struct Energy),comp_erg);
	for(i=0;i<gaslen;i++)
		GasArr[i]=SubEsort[i].PID;//copy the sorted index back

	free(GasE);
	free(SubEsort);
	*P2Len=gaslen;
	*P2PIndex=GasArr;
	
	return gaslen;
}

void make_gassrccat(GASSRCCAT *GSrcCat,GASSUBCAT *GSubCat,GASHALOCAT *GCat,SUBCATALOGUE *SubCat,HBTInt SnapshotNum)
{/*=output:GSrcCat,GSubCat (the same,containing Srccat to be unbinded)
*=input: GSrcCat,GCat,SubCat
* halosrc from current gashalocat
* subhalosrc from its infall gashalocat*/
GASSRCCAT GSrcCatTmp;//Tmp for transfer of SrcCat and then freeing of un-transfered PSubArr;
HBTInt **PSubArr;	// Backup keeps memory of the original PSubArrs for further copy by splitters
HBTInt subid,proid,Npro,Nspl;
HBTInt *sp2pro;
	
	load_sp2pro(SnapshotNum,&Npro,&Nspl,&sp2pro,SUBCAT_DIR);
	if((Npro!=GSrcCat->Nsubs)||(Nspl!=SubCat->Nsplitter))
	{
		fprintf(logfile,"error:pro or spl mismatch:%d=%d,%d=%d\n",Npro,GSrcCat->Nsubs,Nspl,SubCat->Nsplitter);
		exit(1);
	}
	//update GSrcCat
	dump_gassrccat(&GSrcCatTmp,GSrcCat);
	PSubArr=mymalloc(sizeof(HBTInt *)*Npro);
	if(Npro) memcpy(PSubArr,GSrcCatTmp.PSubArr,sizeof(HBTInt *)*Npro);//backup PSubArr
	GSrcCat->Nsubs=SubCat->Nsubs; create_gassrccat(GSrcCat);
	for(subid=0;subid<SubCat->Nsubs-SubCat->NQuasi;subid++)
	{
		if(0==SubCat->SubRank[subid])
		migrate_ghalosrc(GSrcCat,subid,GCat,SubCat->HaloChains[subid].HostID);
		else
		{
		proid=SubCat->HaloChains[subid].ProSubID;
		if(proid>=Npro)
		{
			proid=sp2pro[proid];
			migrate_gsplsrc(GSrcCat,subid,&GSrcCatTmp,proid,PSubArr);//copy from PSubArr 
			//in case the GSrcCatTmp.PSubArr have been or would be transfered to the main-splitters
		}
		else	
		migrate_gsubsrc(GSrcCat,subid,&GSrcCatTmp,proid);
		}
	}
	for(subid=SubCat->Nsubs-SubCat->NQuasi;subid<SubCat->Nsubs;subid++)//quasi halos,no current host
	{
		proid=SubCat->HaloChains[subid].ProSubID;
		if(proid>=Npro)
		{
			proid=sp2pro[proid];
			migrate_gsplsrc(GSrcCat,subid,&GSrcCatTmp,proid,PSubArr);//copy from PSubArr 
			//in case the GSrcCatTmp.PSubArr have been or would be transfered to the main-splitters
		}
		else	
		migrate_gsubsrc(GSrcCat,subid,&GSrcCatTmp,proid);
	}
	myfree(PSubArr);
	erase_gassrccat(&GSrcCatTmp);
	free_sp2pro(sp2pro,Npro,Nspl);

	//make a copy to GSubCat
	GSubCat->Nsubs=SubCat->Nsubs; create_gassubcat(GSubCat);
	memcpy(GSubCat->SubLen,GSrcCat->SubLen,sizeof(HBTInt)*SubCat->Nsubs);
	memcpy(GSubCat->SubOffset,GSrcCat->SubOffset,sizeof(HBTInt)*SubCat->Nsubs);
	for(subid=0;subid<SubCat->Nsubs;subid++)
	{
		GSubCat->PSubArr[subid]=mymalloc(sizeof(HBTInt)*GSubCat->SubLen[subid]);
		memcpy(GSubCat->PSubArr[subid],GSrcCat->PSubArr[subid],sizeof(HBTInt)*GSubCat->SubLen[subid]);
	}
}
void dump_gassrccat(GASSRCCAT *GSrcCatTo,GASSRCCAT *GSrcCatFrom)
{
	memcpy(GSrcCatTo,GSrcCatFrom,sizeof(GASSRCCAT));
	//Nullify SrcCatFrom
	GSrcCatFrom->Nsubs=0;
	GSrcCatFrom->Nids=0;
	GSrcCatFrom->SubLen=NULL;
	GSrcCatFrom->SubOffset=NULL;
	GSrcCatFrom->PSubArr=NULL;
}
void migrate_ghalosrc(GASSRCCAT *GSrcCat, HBTInt subid, GASHALOCAT *Cat,HBTInt hostid)
{
	GSrcCat->SubLen[subid]=Cat->Len[hostid];
	GSrcCat->PSubArr[subid]=mymalloc(sizeof(HBTInt)*Cat->Len[hostid]);
	memcpy(GSrcCat->PSubArr[subid],Cat->PIDorIndex+Cat->Offset[hostid],sizeof(HBTInt)*Cat->Len[hostid]);
}
void migrate_gsubsrc(GASSRCCAT *GSrcCatTo, HBTInt subid, GASSRCCAT *GSrcCatFrom,HBTInt proid)
{
	GSrcCatTo->SubLen[subid]=GSrcCatFrom->SubLen[proid];
	GSrcCatTo->PSubArr[subid]=GSrcCatFrom->PSubArr[proid];GSrcCatFrom->PSubArr[proid]=NULL;
}
void migrate_gsplsrc(GASSRCCAT *GSrcCatTo, HBTInt subid, GASSRCCAT *GSrcCatFrom,HBTInt proid,HBTInt **PSubArr)
{
	GSrcCatTo->SubLen[subid]=GSrcCatFrom->SubLen[proid];
	GSrcCatTo->PSubArr[subid]=mymalloc(sizeof(HBTInt)*GSrcCatFrom->SubLen[proid]);
	memcpy(GSrcCatTo->PSubArr[subid],PSubArr[proid],sizeof(HBTInt)*GSrcCatFrom->SubLen[proid]);
}

#define NDIV_gas 200
static HBTInt hoc_gas[NDIV_gas][NDIV_gas][NDIV_gas],ll_gas[NP_GAS];
static HBTReal range_gas[3][2], step_gas[3];
void makell_gas()
{
	HBTReal (*pos)[3];
	HBTInt i,j,grid[3],np;
	
	pos=Gdat.Pos;
	np=NP_GAS;
	printf("creating linked list..\n");
	
	/*determining enclosing cube*/
	for(i=0;i<3;i++)
		for(j=0;j<2;j++)
			range_gas[i][j]=pos[0][i];
	for(i=1;i<np;i++)
		for(j=0;j<3;j++)
		{
			if(pos[i][j]<range_gas[j][0])
				range_gas[j][0]=pos[i][j];
			else if(pos[i][j]>range_gas[j][1])
				range_gas[j][1]=pos[i][j];
		}
	for(j=0;j<3;j++)
		step_gas[j]=(range_gas[j][1]-range_gas[j][0])/NDIV_gas;
	
	/*initialize hoc*/
	HBTInt *phoc=&(hoc_gas[0][0][0]);
	for(i=0;i<NDIV_gas*NDIV_gas*NDIV_gas;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((pos[i][j]-range_gas[j][0])/step_gas[j]);
			if(grid[j]<0) 
				grid[j]=0;
			else if(grid[j]>=NDIV_gas)
				grid[j]=NDIV_gas-1;
		}			
		ll_gas[i]=hoc_gas[grid[0]][grid[1]][grid[2]];
		hoc_gas[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
}

void collect_gas_particles(HBTInt subid,SUBCATALOGUE*SubCat,HBTInt *P2GasSrcLen,HBTInt **P2GasSrc)
{
	HBTInt gaslen,i,j,k,pid,subbox_grid[3][2],maxlen,*GasSrc;
	HBTReal *cen,rvir,rscale,dr;
	//~ extern HBTReal range[3][2],step[3],hoc[][][],ll[],Rtidal[];
	
	gaslen=0;
	if(SubCat->SubLen[subid])
	{
	maxlen=SubCat->SubLen[subid]*2;
	maxlen=((maxlen<NP_GAS)?maxlen:NP_GAS);
	GasSrc=mymalloc(sizeof(HBTInt)*maxlen);
	cen=SubCat->Property[subid].CoM;
	rvir=comoving_virial_radius(SubCat->SubLen[subid]);
	rscale=rvir*GasCollect_ScaleRelax;
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rscale-range_gas[i][0])/step_gas[i]);
	if(subbox_grid[i][0]<0)subbox_grid[i][0]=0;
	subbox_grid[i][1]=floor((cen[i]+rscale-range_gas[i][0])/step_gas[i]);
	if(subbox_grid[i][1]>=NDIV_gas)subbox_grid[i][1]=NDIV_gas-1;
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=hoc_gas[i][j][k];
				while(pid>=0)
				{
					dr=distance(Gdat.Pos[pid],cen);
					if(dr<rscale)
					{
						GasSrc[gaslen]=pid;
						gaslen++;
						if(gaslen>=maxlen)
						{
							maxlen*=2;
							GasSrc=realloc(GasSrc,sizeof(HBTInt)*maxlen);
						}
					}
					pid=ll_gas[pid];
				}
			}
			GasSrc=realloc(GasSrc,sizeof(HBTInt)*gaslen);
	}
	else
	{
		GasSrc=NULL;
	}
	*P2GasSrcLen=gaslen;
	*P2GasSrc=GasSrc;
}
