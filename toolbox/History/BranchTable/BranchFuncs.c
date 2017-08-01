#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "BranchFuncs.h"
HBTReal * load_Vmax(HBTInt Nsnap)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	HBTReal *vmax;
        
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
        fseek(fp, sizeof(HBTReal)*Nsubs, SEEK_CUR);
// 	fread(rmax,sizeof(HBTReal),Nsubs,fp);
        vmax=mymalloc(sizeof(HBTReal)*Nsubs);
	fread(vmax,sizeof(HBTReal),Nsubs,fp);
	fseek(fp,sizeof(HBTReal)*Nsubs*5,SEEK_CUR);
//         fread(rhalf,sizeof(HBTReal),Nsubs,fp);
	//~ fread(rsig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r2sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r3sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(rpoisson,sizeof(HBTReal),Nsubs,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
	return vmax;
}
void section_ID2Index(CrossSection *sec)
{
  HBTInt NodeID;
  if(0==PIDHash.np)
  {
	  fprintf(logfile,"call fill_PIDHash() before section_ID2Index()!\n");
	  exit(1);
  }
  #pragma omp parallel for
  for(NodeID=0;NodeID<sec->NumNode;NodeID++)
	sec->Node[NodeID].MstBndID=lookup_ID2Ind(sec->Node[NodeID].MstBndID);
}
void section_Ind2ID(CrossSection *sec)
{
  HBTInt NodeID;
  #pragma omp parallel for
  for(NodeID=0;NodeID<sec->NumNode;NodeID++)
	sec->Node[NodeID].MstBndID=lookup_Ind2ID(sec->Node[NodeID].MstBndID);
}
void section_init(CrossSection * sec)
{
  sec->NumNode=0;
  sec->NumNodeAlloc=NumNodeAllocMin;
  sec->Node=mymalloc(sizeof(BranchNode)*sec->NumNodeAlloc);
}
void section_free(CrossSection *sec)
{
  sec->NumNode=0;
  sec->NumNodeAlloc=0;
  myfree(sec->Node);
}
void save_cross_section(HBTInt Nsnap, CrossSection *sec, char *outdir)
{
  FILE *fp;
  char buf[1024];//outdir[1024];
//   sprintf(outdir, "%s/BranchTable",SUBCAT_DIR);
  mkdir(outdir,0755);
  sprintf(buf, "%s/CrossSection_%03d", outdir, (int)Nsnap);
  myfopen(fp,buf,"w");
  fwrite(&sec->NumNode, sizeof(HBTInt), 1, fp);
  fwrite(sec->Node, sizeof(BranchNode), sec->NumNode, fp);
  fwrite(&sec->NumNode, sizeof(HBTInt), 1, fp);
  fclose(fp);
}

void load_cross_section(HBTInt Nsnap, CrossSection * sec, char *path)
{
  if(Nsnap<IniSnap)//load nothing
  {
	section_init(sec);
	return;
  }
  
  FILE *fp;
  char buf[1024];
  sprintf(buf, "%s/CrossSection_%03d", path, (int)Nsnap);
  myfopen(fp,buf,"r");
  fread(&sec->NumNode, sizeof(HBTInt), 1, fp);
  sec->NumNodeAlloc=(sec->NumNode<NumNodeAllocMin?NumNodeAllocMin:sec->NumNode);
  sec->Node=mymalloc(sizeof(BranchNode)*sec->NumNodeAlloc);
  fseek(fp, sizeof(BranchNode)*sec->NumNode, SEEK_CUR);
  HBTInt n;
  fread(&n, sizeof(HBTInt), 1, fp);
  fclose(fp);
  if(n!=sec->NumNode)
  {
	fprintf(logfile, "Error: numbers of nodes do not match in %s\n"HBTIFMT","HBTIFMT"\nFile corruption?\n", buf, sec->NumNode, n);
	exit(1);
  }
  
#define NUM_THREAD 10
#pragma omp parallel num_threads(NUM_THREAD) 
  {
    int ithread=omp_get_thread_num();
    HBTInt nread=sec->NumNode/NUM_THREAD;
    HBTInt offset=nread*ithread;
    if(ithread==NUM_THREAD-1)
        nread=sec->NumNode-offset;
    FILE *fpthis;
    myfopen(fpthis, buf, "r");
    fseek(fpthis, sizeof(HBTInt)+sizeof(BranchNode)*offset, SEEK_SET);
    fread(sec->Node+offset, sizeof(BranchNode), nread, fpthis);
    fclose(fpthis);
  }

}
void get_phys_vel(HBTxyz vel, HBTInt PID)
{
  #ifdef VEL_INPUT_PHYSICAL
  memcpy(vel, Pdat.Vel+PID, sizeof(HBTxyz));
  #else
  int i;
  for(i=0;i<3;i++)
	vel[i]=Pdat.Vel[PID][i]*sqrt(header.time);
  #endif
}
void section_fill_new_node(CrossSection *sec, HBTInt NodeID, HBTInt SubID, SUBCATALOGUE *SubCat, HBTReal *Vmax)
{
  if(NodeID>=sec->NumNodeAlloc)//the section should be large enough to hold the new node
  {
	fprintf(logfile, "Error: NodeID=%ld over flow (%ld). Reallocate CrossSection.Node before proceeding.\n", (long)NodeID, (long)sec->NumNodeAlloc);
	exit(1);
  }
  BranchNode *node=sec->Node+NodeID;
  node->BranchID=NodeID;
  node->SubID=SubID;
  node->HostID=SubCat->HaloChains[SubID].HostID;
  node->SubRank=SubCat->SubRank[SubID];
  if(SubCat->SubLen[SubID])
  {
	node->MstBndID=SubCat->PSubArr[SubID][0];
	memcpy(node->MstBndPos, Pdat.Pos+node->MstBndID, sizeof(HBTxyz));
	get_phys_vel(node->MstBndVel, node->MstBndID);
  }
  else
  {
	node->MstBndID=-1;
	int i;
	for(i=0;i<3;i++)
	{
	  node->MstBndPos[i]=0.;
	  node->MstBndVel[i]=0.;
	}
  }
  node->NpBnd=SubCat->SubLen[SubID];
  node->NpBndPeak=node->NpBnd;
  node->SnapNumPeak=Pdat.Nsnap;
  node->Vmax=Vmax[SubID];
  node->VmaxPeak=node->Vmax;
  node->SnapNumVpeak=Pdat.Nsnap;
}
void section_create_node(CrossSection * sec, HBTInt NodeID, HBTInt SubID)
{
  if(NodeID==sec->NumNodeAlloc)
  {
	sec->NumNodeAlloc*=2;
	sec->Node=realloc(sec->Node, sizeof(BranchNode)*sec->NumNodeAlloc);
  }
  sec->Node[NodeID].BranchID=NodeID;
  sec->Node[NodeID].SubID=SubID;
}
void node_fill(BranchNode *node, SUBCATALOGUE *SubCat, HBTReal *Vmax)
{
  HBTInt SubID=node->SubID;
  node->HostID=SubCat->HaloChains[SubID].HostID;
  node->SubRank=SubCat->SubRank[SubID];
  if(SubCat->SubLen[SubID])
  {
	node->MstBndID=SubCat->PSubArr[SubID][0];
	memcpy(node->MstBndPos, Pdat.Pos+node->MstBndID, sizeof(HBTxyz));
	get_phys_vel(node->MstBndVel, node->MstBndID);
  }
  else
  {
	node->MstBndID=-1;
	int i;
	for(i=0;i<3;i++)
	{
	  node->MstBndPos[i]=0.;
	  node->MstBndVel[i]=0.;
	}
  }
  node->NpBnd=SubCat->SubLen[SubID];
  node->NpBndPeak=node->NpBnd;
  node->SnapNumPeak=Pdat.Nsnap;
  node->Vmax=Vmax[SubID];
  node->VmaxPeak=node->Vmax;
  node->SnapNumVpeak=Pdat.Nsnap;
}
void node_update(BranchNode *node, HBTInt *pro2dest, SUBCATALOGUE *SubCat, CATALOGUE *Cat, HBTReal *Vmax)
{
  HBTInt DestID;
  DestID=pro2dest[node->SubID]; //this maps -1 to -1.
  if(DestID<0)//disrupted
  {
	//BranchID,MstBndID,NpBndPeak,SnapBndPeak does not change.
	node->SubID=-1;
	node->SubRank=-1; //disrupted subhalo
	if(node->MstBndID<0) //no bound particle
	  node->HostID=-1;
	else
	{
	  node->HostID=Cat->ID2Halo[node->MstBndID];
	  memcpy(node->MstBndPos, Pdat.Pos+node->MstBndID, sizeof(HBTxyz));
	  get_phys_vel(node->MstBndVel, node->MstBndID);
	}
	node->NpBnd=0;
        node->Vmax=0;
  }
  else
  {
	//BranchID does not change
	node->SubID=DestID;
	node->SubRank=SubCat->SubRank[DestID];//there could be cases when NpBnd=0 but SubRank, SubID and MstBnd info still make sense
	node->HostID=SubCat->HaloChains[DestID].HostID;
	if(SubCat->SubLen[DestID])//only update mstbndID if the dest has non-zero bound mass
	  node->MstBndID=SubCat->PSubArr[DestID][0];//PSubArr already converted to index
	memcpy(node->MstBndPos, Pdat.Pos+node->MstBndID, sizeof(HBTxyz));
	get_phys_vel(node->MstBndVel, node->MstBndID);
	node->NpBnd=SubCat->SubLen[DestID];
	if(node->NpBndPeak<node->NpBnd)
	{
	  node->NpBndPeak=node->NpBnd;
	  node->SnapNumPeak=Pdat.Nsnap;
	}
	node->Vmax=Vmax[DestID];
        if(node->VmaxPeak<node->Vmax)
        {
            node->VmaxPeak=node->Vmax;
            node->SnapNumVpeak=Pdat.Nsnap;
        }
  }
}
void parse_snap_args(HBTInt SnapRange[2], int argc, char **argv)
{
	if(argc==1)
	{
		SnapRange[0]=IniSnap;
		SnapRange[1]=MaxSnap-1;
	}
	 else if(argc==2)
	 {
		 SnapRange[0]=atoi(argv[1]);
		 SnapRange[1]=atoi(argv[1]);
	 }
	 else if(argc==3)
	 {
		 SnapRange[0]=atoi(argv[1]);
		 SnapRange[1]=atoi(argv[2]);
	 }
	 if((argc>3)||SnapRange[0]<IniSnap||SnapRange[0]>=MaxSnap||SnapRange[1]<IniSnap||SnapRange[1]>=MaxSnap||SnapRange[1]<SnapRange[0])
	  {
		 printf("Usage: %s [SnapShot_begin [SnapShot_end]]\n\
				\tSnapShot_begin and SnapShot_end being integers in the range IniSnap~MaxSnap-1,specifying the snap range to be run\n\
				\t if no param is specified, run over all the snapshots\n\
				\t if only SnapShot_begin is specified, Snapshot_end is set to SnapShot_begin\n",argv[0]);
		 exit(1);
	 }
}