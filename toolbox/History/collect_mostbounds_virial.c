/*collect all the mostbound particles inside 2*Rv
 * tag their infall according to: 1) Mmax 2) Rvir cross (first and last)
 * tag direct/indirect infall according to whether Mmax is directly hosted by the central
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "hdf_util.h"

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "history_vars.h"
#include "history_proto.h"

#define VIRTYPE 1

struct HistoryShards HistoryRevShard;
HBTInt * (MainSubID[MaxSnap]);
EVOLUTIONCAT EvoCatPre,EvoCat;
HALOSIZE *(halosize[MaxSnap]);
HBTInt * (Sub2Hist[MaxSnap]);
HBTxyz * (SubCoM[MaxSnap]);

struct PList
{
  HBTInt CenSubID;
  HBTReal Cen[3], VCen[3], Rmax;
  
  HBTInt np;
  HBTInt *PIndex;
  float (* pos)[3],(* vel)[3];
  int * mass, *IsDirectInfall, (*snapTVV)[3];
  int (*massTVV)[3];
};

MBDCATALOGUE MbdCat;
SUBCATALOGUE SubCat;
extern void load_evocat_raw(EVOLUTIONCAT *EvoCat);
extern void collect_particles(struct PList * p);
extern void fill_particle_nodes(struct PList *p);
extern void dump_particles_hdf(char *outfile, struct PList *p);
int main(int argc, char **argv)
{
  
  HBTInt Nsnap,i;
  logfile=stdout;
  
  #pragma omp parallel for 
  for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
  {
  HBTInt Ngroups, Nsubs;
  load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
  Ngroups=read_Ngroups(Nsnap);
  halosize[Nsnap]=mymalloc(sizeof(HALOSIZE)*Ngroups);
  load_halo_size(halosize[Nsnap],Ngroups,Nsnap);
  load_mainsubid(Nsnap,MainSubID+Nsnap);
  read_subpos(Nsnap,SubCoM+Nsnap);
  printf(HBTIFMT".",Nsnap);fflush(stdout);
  }
  load_particle_header(MaxSnap-1, SNAPSHOT_DIR);

  load_evocat_raw(&EvoCatPre);
  load_evocat_rev(&EvoCat,SUBCAT_DIR);
  load_historyshards(&HistoryRevShard);
  
  Nsnap=MaxSnap-1;
  load_mbd_catalogue(Nsnap, &MbdCat);
  printf("%d-%d-%d particles\n", MbdCat.NSubs+MbdCat.NOrphans, MbdCat.GrpLen_Sub[0], MbdCat.GrpLen_Orphan[0]);
  load_sub_table(Nsnap, &SubCat, SUBCAT_DIR);
  printf("%d-%d subs\n", SubCat.Nsubs, SubCat.GrpLen_Sub[0]);
  
  struct PList data;
  data.CenSubID=0;
  data.Rmax=3.*halosize[Nsnap][0].Rvir[VIRTYPE];
  for(i=0;i<3;i++)
  {//or MbdCat.Node[CenSubID].Pos;
	data.Cen[i]=SubCat.Property[data.CenSubID].CoM[i];
	data.VCen[i]=SubCat.Property[data.CenSubID].VCoM[i];
  }
  collect_particles(&data);
  fill_particle_nodes(&data);
  char outfile[1024];
  sprintf(outfile, "%s/anal/MbdInfall_VIR%d.hdf5", SUBCAT_DIR, VIRTYPE);
  dump_particles_hdf(outfile, &data);
  
  return 0;
}

void load_evocat_raw(EVOLUTIONCAT *EvoCat)
{
	HBTInt HistID,Nsnap;
	SubNode *Member;
	load_history_pre(EvoCat,SUBCAT_DIR);
	EvoCat->PartMass=header.mass[1];;
	for(HistID=0;HistID<EvoCat->NHist;HistID++)
	{
		for(Nsnap=EvoCat->History[HistID].SnapBirth;Nsnap<EvoCat->History[HistID].SnapDeath;Nsnap++)
		{
			Member=EvoCat->History[HistID].Member-EvoCat->History[HistID].SnapBirth+Nsnap;
			HBTInt HostID;
			HostID=Member->HostID;
			if(HostID<0)//quasi halo
			{
			Member->Mhost=0;
			//~ Member->Chost=0;
			}
			else
			{
			//~ Member->Mhost=Cat.Len[HostID];
			Member->Mhost=halosize[Nsnap][HostID].Mvir[0];
			//~ Member->Chost=halocon[Nsnap][HostID];
			}
		}
	}
}
HBTInt get_mainsubid(HBTInt Nsnap,HBTInt hostid)
{
	 if(hostid<0)
  {
	  printf("error: wrong hostid to read for mainsub\n");
	  exit(1);
  }
  return MainSubID[Nsnap][hostid];
}
HBTInt is_direct_satellite(HBTInt CenHistID,HBTInt HistID, HBTInt Nsnap)//return 1 if the current host is grpmainid
{
	HBTInt HostID, CenID;

	if(HistID==CenHistID)	
	  return 0;//do not count the main-branch itself
	
	if(0==HistoryRevShard.NBirth[HistID]) 
	  return -1; //never infalled, always isolate
	
	HostID=GetMember(&EvoCat,HistID,Nsnap)->HostID;	
// 	if(HostID<0) return -1; //quasi halo and never infall
	CenID=GetMember(&EvoCat,CenHistID,Nsnap)->HostID;
	if(CenID<0||CenID!=HostID) return 0;
	
	return 1;
}
HBTReal * GetMbdPos(HBTInt subid,void *data)
{//pass SubCat.Property to data
	return ((MBDNode *) data)[subid].Pos;
}
void collect_particles(struct PList * p)
{
  LINKLIST ll;
  make_linklist(&ll, MbdCat.NSubs+MbdCat.NOrphans, 50, MbdCat.Nodes, GetMbdPos, 0);
  p->np=MbdCat.NSubs+MbdCat.NOrphans;
  p->PIndex=linklist_search_sphere(&ll, p->Rmax, p->Cen, &p->np);
  printf("found %d particles\n", (int)p->np);
  free_linklist(&ll);
}
HBTInt get_Mmax_Node(HBTInt HistID)
{
  HBTInt i,imax=0,Mmax=0,m;
  for(i=0;i<EvoCat.HistLen[HistID];i++)
  {
	m=EvoCat.History[HistID].Member[i].Mdm;
	if(m>Mmax)
	{
	  imax=i;
	  Mmax=m;
	}
  }
  return imax;
}
void get_infall_snap(HBTInt CenHistID, HBTInt HistID, int SnapInfall[2])
{
  HBTInt Nsnap, SnapBirth, InHalo=0;
  SnapInfall[0]=-1;
  SnapInfall[1]=-1;
#define MAX(a,b) ((a)>(b)?(a):(b))
  SnapBirth=MAX(EvoCat.History[HistID].SnapBirth, EvoCat.History[CenHistID].SnapBirth);
  for(Nsnap=SnapBirth;Nsnap<EvoCat.History[HistID].SnapDeath;Nsnap++)
  {
	SubNode *CenNode=GetMember(&EvoCat, CenHistID, Nsnap);
	SubNode *SatNode=GetMember(&EvoCat, HistID, Nsnap);
	HBTReal r;
	if(CenNode!=NULL&&SatNode!=NULL)
	{
	  HALOSIZE * halo=&(halosize[Nsnap][CenNode->HostID]);
	  if(0==halo->flag_badvir[VIRTYPE]&&0==halo->flag_fakehalo)//skip bad halos
	  {
		r=distance(SubCoM[Nsnap][CenNode->SubID], SubCoM[Nsnap][SatNode->SubID]);
		if(r<halo->Rvir[VIRTYPE])
		{
		  if(SnapInfall[0]<0) SnapInfall[0]=Nsnap;//first infall
		  if(!InHalo) //just falled in.
		  {
			SnapInfall[1]=Nsnap;// record for last infall
			InHalo=1;
		  }
		}
		else//outside
		{
		  if(InHalo)  InHalo=0;//toggle
		}
	  }
	}
  }
}
void fill_particle_nodes(struct PList *p)
{
  int i,np=p->np;
  p->pos=mymalloc(sizeof(float)*3*p->np);
  p->vel=mymalloc(sizeof(float)*3*p->np);
  p->mass=mymalloc(sizeof(int)*p->np);
  p->IsDirectInfall=mymalloc(sizeof(int)*np);//level-0 satellite at Mmax
  p->snapTVV=mymalloc(sizeof(float)*np*3); //Mmax, FirstInfall, LastInfall
  p->massTVV=mymalloc(sizeof(int)*np*3);
  
  HBTInt CenHistID=Sub2Hist[MaxSnap-1][p->CenSubID];
#pragma omp parallel for
  for(i=0;i<np;i++)
  {
	  MBDNode *Node;
	  HBTInt j,NodeID, HistID, imax, ShardID;
	  NodeID=p->PIndex[i];
	  Node=MbdCat.Nodes+NodeID;
	  for(j=0;j<3;j++)
	  {
		p->pos[i][j]=Node->Pos[j]-p->Cen[j];
		p->vel[i][j]=Node->Vel[j]-p->VCen[j];
	  }
	  if(NodeID<SubCat.Nsubs)
		p->mass[i]=SubCat.SubLen[NodeID];
	  else
		p->mass[i]=0;
	  
	  HistID=Node->HistID;
	  imax=get_Mmax_Node(HistID);
	  p->massTVV[i][0]=EvoCat.History[HistID].Member[imax].Mdm;
	  p->snapTVV[i][0]=EvoCat.History[HistID].SnapBirth+imax;
	  p->IsDirectInfall[i]=is_direct_satellite(CenHistID, HistID, p->snapTVV[i][0]);
	  get_infall_snap(CenHistID, HistID, p->snapTVV[i]+1);
	  for(j=1;j<3;j++)
	  {
		if(p->snapTVV[i][j]>=0)
		  p->massTVV[i][j]=GetMember(&EvoCat, HistID, p->snapTVV[i][j])->Mdm;
	  }
  }
}

void dump_particles_hdf(char *outfile, struct PList *p)
{
    hid_t    file_id,group_id;
    herr_t      status;
    size_t     i, j,nread=0;
    hsize_t dims[2];
	int np=p->np;
	
    file_id = H5Fcreate (outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file
	
    dims[0]=1;
	float pmass=header.mass[1];
    status = H5LTmake_dataset(file_id,"/PartMass", 1, dims, H5T_NATIVE_FLOAT, &pmass); 
    dims[0]=np;
    dims[1]=3;
    status = H5LTmake_dataset(file_id,"/x",2,dims,H5T_NATIVE_FLOAT,p->pos);
    status = H5LTmake_dataset(file_id,"/v",2,dims,H5T_NATIVE_FLOAT,p->vel);
	float x0[3],v0[3];
	for(j=0;j<3;j++)
	{
	  x0[j]=p->Cen[j];
	  v0[j]=p->VCen[j];
	}
	float rmax=p->Rmax;
	H5LTset_attribute_float(file_id, "/x", "x0", x0, 3);
	H5LTset_attribute_float(file_id, "/v", "v0", v0, 3);
	H5LTset_attribute_float(file_id, "/x", "rmax", &rmax, 1);
	
	dims[0]=np;
    dims[1]=3;
	status = H5LTmake_dataset(file_id,"/DirectInfall",1,&dims[0],H5T_NATIVE_INT,p->IsDirectInfall);//0, 1, -1(neverinfalled) for Mmax.
    status = H5LTmake_dataset(file_id,"/mass",1,&dims[0],H5T_NATIVE_INT,p->mass);//in particles
    status = H5LTmake_dataset(file_id,"/snapTVV",2,dims,H5T_NATIVE_INT,p->snapTVV);//SnapTidal, SnapFirstRvInfall, SnapLastRvInfall
	status = H5LTmake_dataset(file_id,"/massTVV",2,dims,H5T_NATIVE_INT,p->massTVV);//in particles
		
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[0]: (%g,%g,%g)\n",p->pos[0][0],p->pos[0][1],p->pos[0][2]);
}