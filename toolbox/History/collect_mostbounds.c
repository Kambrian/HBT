/*collect all the mostbound particles inside 2*Rv
 * for each particle, list their pos, vel, current mass, SnapInfall(fof), SnapRvir, SnapTidal, IsDirectInfall, mass (bound and Mvir?) at these times, orbital parameters at these times
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

struct HistoryShards HistoryRevShard;
HBTInt * (MainSubID[MaxSnap]);
EVOLUTIONCAT EvoCatPre,EvoCat;
HALOSIZE *(halosize[MaxSnap]);
HBTInt * (Sub2Hist[MaxSnap]);

struct PList
{
  HBTInt CenSubID;
  HBTReal Cen[3], VCen[3], Rmax;
  
  HBTInt np;
  HBTInt *PIndex;
  float (* pos)[3],(* vel)[3];
  int * mass, *IsDirectInfall, (*snapTV)[2] ;
  float (*massTV)[2], (*kappaTV)[2], (*j2TV)[2];
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
  data.Rmax=3.*halosize[Nsnap][0].Rvir[0];
  for(i=0;i<3;i++)
  {//or MbdCat.Node[CenSubID].Pos;
	data.Cen[i]=SubCat.Property[data.CenSubID].CoM[i];
	data.VCen[i]=SubCat.Property[data.CenSubID].VCoM[i];
  }
  collect_particles(&data);
  fill_particle_nodes(&data);
  char outfile[1024];
  sprintf(outfile, "%s/anal/MbdInfall.hdf5", SUBCAT_DIR);
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
HBTInt is_direct_grp_family(HBTInt grpmainid,HBTInt HistID,HBTInt ShardID)//return 1 if this shard falls directly to the central halo
{
	HBTInt SnapRvir,HostID;

	if(HistID==Sub2Hist[MaxSnap-1][grpmainid])	return -1;//do not count the main-branch itself

	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	if(SnapRvir<0)	return -1; //no host; never infalled
	
	HostID=GetMember(&EvoCat,HistID,SnapRvir)->HostID;
	HistID=Sub2Hist[SnapRvir][get_mainsubid(SnapRvir,HostID)];
	if(HistID==Sub2Hist[MaxSnap-1][grpmainid])
	return 1;
	else
	return 0;
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
void fill_particle_nodes(struct PList *p)
{
	int i,j,np=p->np;
	
    p->pos=mymalloc(sizeof(float)*3*np);
    p->vel=mymalloc(sizeof(float)*3*np);
	p->mass=mymalloc(sizeof(int)*np);
	p->IsDirectInfall=mymalloc(sizeof(int)*np);
	p->snapTV=mymalloc(sizeof(int)*np*2);
	p->massTV=mymalloc(sizeof(float)*np*2);
	p->kappaTV=mymalloc(sizeof(float)*np*2);
	p->j2TV=mymalloc(sizeof(float)*np*2);
	
    for(i=0;i<np;i++)
    {
	  MBDNode *Node;
	  int NodeID, HistID, ShardID;
	  NodeID=p->PIndex[i];
	  Node=MbdCat.Nodes+NodeID;
	  HistID=Node->HistID;
	  ShardID=HistoryRevShard.NBirth[HistID]-1;//the final shard
      
	  for(j=0;j<3;j++)
      {
		p->pos[i][j]=Node->Pos[j]-p->Cen[j];
		p->vel[i][j]=Node->Vel[j]-p->VCen[j];
      }
      if(NodeID<SubCat.Nsubs)
		p->mass[i]=SubCat.SubLen[NodeID];
	  else
		p->mass[i]=0;
	  
	  if(ShardID>=0)
	  {
		struct ShardParam *Par;
		Par=HistoryRevShard.Par[HistID]+ShardID;
		p->snapTV[i][0]=Par->SnapTidal;
		p->snapTV[i][1]=Par->SnapRvir;
		for(j=0;j<2;j++)
		{
		  p->massTV[i][j]=Par->Mhost[j]*Par->Mrate[j];
		  p->kappaTV[i][j]=Par->Kappa[j];
		  p->j2TV[i][j]=Par->j2[j];
		}
		p->IsDirectInfall[i]=is_direct_grp_family(p->CenSubID, HistID, ShardID);
	  }
	  else //never infalled
	  {
		p->snapTV[i][0]=-1;
		p->snapTV[i][1]=-1;
		for(j=0;j<2;j++)
		{
		  p->massTV[i][j]=p->mass[i];//current mass
		  p->kappaTV[i][j]=0;
		  p->j2TV[i][j]=0;
		}
		p->IsDirectInfall[i]=-1;//never infalled
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
    dims[1]=2;
	status = H5LTmake_dataset(file_id,"/DirectInfall",1,&dims[0],H5T_NATIVE_INT,p->IsDirectInfall);//-1(never infalled), 0, 1.
    status = H5LTmake_dataset(file_id,"/mass",1,&dims[0],H5T_NATIVE_INT,p->mass);//in particles
    status = H5LTmake_dataset(file_id,"/snapTV",2,dims,H5T_NATIVE_INT,p->snapTV);//SnapTidal, SnapRvir
	status = H5LTmake_dataset(file_id,"/massTV",2,dims,H5T_NATIVE_FLOAT,p->massTV);//in particles
	status = H5LTmake_dataset(file_id,"/kappaTV",2,dims,H5T_NATIVE_FLOAT,p->kappaTV);
	status = H5LTmake_dataset(file_id,"/j2TV",2,dims,H5T_NATIVE_FLOAT,p->j2TV);
		
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[0]: (%g,%g,%g)\n",p->pos[0][0],p->pos[0][1],p->pos[0][2]);
}