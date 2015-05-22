/*collect all the mostbound particles inside 2*Rv, from DTrees on SubFind
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

#ifdef CONVERT_LENGTH_MPC_KPC
#define LSCALE 1e-3
#else 
#define LSCALE 1.
#endif

struct PList
{
  HBTInt CenBranchID;
  HBTReal Cen[3], VCen[3], Rmax;
  
  HBTInt np;
  HBTInt *PIndex;
  float (* pos)[3],(* vel)[3];
  int * mass, *IsDirectInfall, (*snapTVV)[3];
  int (*massTVV)[3];
};
struct TreeNode
{
  int SnapNum;
  int IsFoFCentre;
  HBTInt MBDID;
  HBTInt nP;
  HBTInt nP_infall;
  HBTxyz Pos;
};
struct TreeBranch
{
  HBTInt NNode;
  struct TreeNode *Node;
};
struct Tree
{
  HBTInt NBranch;
  struct TreeBranch *Branch;
};
struct Tree Tree;
extern void load_tree(struct Tree *Tree, char *infile);
extern void free_tree(struct Tree *Tree);
extern int get_mostmassive_branch(struct Tree *Tree);
extern HBTInt get_Mmax_Node(HBTInt BranchID);
extern void get_infall_nodes(HBTInt BranchID, HBTInt InfallNodeID[2]);
extern void collect_particles(struct PList * p, HBTInt *MBDList, HBTInt np);
extern void fill_particle_nodes(struct PList *p, HBTInt *MBDList);
extern void dump_particles_hdf(char *outfile, struct PList *p);
int main(int argc, char **argv)
{
  
  HBTInt Nsnap,i;
  logfile=stdout;
  
  struct PList data;
//   load_tree(&Tree, "/gpfs/data/jch/Aquarius/Trees/Aq-A/1/trees/treedir_060/tree_060.0.hdf5");
//   load_tree(&Tree, "/gpfs/data/jch/Aquarius/Trees/Aq-A/2/trees/treedir_127/tree_127.0.hdf5");
//   load_tree(&Tree, "/gpfs/data/jch/Aquarius/Trees/Aq-A/3/trees/treedir_063/tree_063.0.hdf5");
//   load_tree(&Tree, "/gpfs/data/jch/Aquarius/Trees/Aq-A/4/trees/treedir_127/tree_127.0.hdf5");
//   load_tree(&Tree, "/gpfs/data/jch/Aquarius/Trees/Aq-A/5/trees/treedir_127/tree_127.0.hdf5");
//   load_tree(&Tree, "/gpfs/data/jch/Aquarius/Trees/Aq-B/4/trees/treedir_127/tree_127.0.hdf5");
load_tree(&Tree, "/gpfs/data/jch/Aquarius/Trees/Aq-F/2/trees/treedir_097/tree_097.0.hdf5");
// #define RVIR 179.4e-3 //Mpc/h, A
// #define RVIR 137.86e-3 //B
#define RVIR 180e-3 //Mpc/h, conservative.
//   load_tree(&Tree, "/gpfs/data/jch/Phoenix/Trees/pha-2/trees/treedir_071/tree_071.0.hdf5");
// #define RVIR 1.414 //PhA
// #define RVIR 1.526 //B
// #define RVIR 1.332 //C
// #define RVIR 1.386 //D
// #define RVIR 1.369 //E
// #define RVIR 1.509 //F
// #define RVIR 1.704 //G
// #define RVIR 2.411 //I
  data.CenBranchID=get_mostmassive_branch(&Tree);
  
  Nsnap=MaxSnap-1;
  load_particle_data(Nsnap, SNAPSHOT_DIR);
  
  HBTInt *MBDList=malloc(sizeof(HBTInt)*Tree.NBranch);
  for(i=0;i<Tree.NBranch;i++)
	MBDList[i]=Tree.Branch[i].Node[0].MBDID;
  fill_PIDHash();
  fresh_ID2Index(MBDList, Tree.NBranch);
  free_PIDHash();
  
  data.Rmax=3.*RVIR;
  for(i=0;i<3;i++)
  {
	data.Cen[i]=Pdat.Pos[MBDList[data.CenBranchID]][i]*LSCALE;
	data.VCen[i]=Pdat.Vel[MBDList[data.CenBranchID]][i];
  }
  collect_particles(&data, MBDList, Tree.NBranch);
  
  fill_particle_nodes(&data, MBDList);
  char outfile[1024];
  sprintf(outfile, "%s/anal/MbdInfall.DTree.hdf5", SUBCAT_DIR);
  dump_particles_hdf(outfile, &data);
  
  return 0;
}

HBTReal * GetMbdPos(HBTInt i, void *data)
{
	return Pdat.Pos[((HBTInt *) data)[i]];
}
void collect_particles(struct PList * p, HBTInt *MBDList, HBTInt np)
{//return the list of branches, not MBDIndex.
  LINKLIST ll;
  make_linklist(&ll, np, 50, MBDList, GetMbdPos, 0);
  p->np=np;
  HBTReal Cen[3];
  HBTInt i;
  for(i=0;i<3;i++) Cen[i]=p->Cen[i]/LSCALE; //convert to raw unit of Pdat
  p->PIndex=linklist_search_sphere(&ll, p->Rmax/LSCALE, Cen, &p->np);
  printf("found %d particles\n", (int)p->np);
  free_linklist(&ll);
}
void fill_particle_nodes(struct PList *p, HBTInt *MBDList)
{
  int i,np=p->np;
  p->pos=mymalloc(sizeof(float)*3*p->np);
  p->vel=mymalloc(sizeof(float)*3*p->np);
  p->mass=mymalloc(sizeof(int)*p->np);
  p->IsDirectInfall=mymalloc(sizeof(int)*np);//level-0 satellite at Mmax
  p->snapTVV=mymalloc(sizeof(float)*np*3); //Mmax, FirstInfall, LastInfall
  p->massTVV=mymalloc(sizeof(int)*np*3);

HBTInt FinalSnap=Tree.Branch[p->CenBranchID].Node[0].SnapNum;
#pragma omp parallel for
  for(i=0;i<np;i++)
  {
	int j,k;
	HBTInt BranchID=p->PIndex[i];
	HBTInt pid=MBDList[BranchID];
	  for(j=0;j<3;j++)
	  {
		p->pos[i][j]=Pdat.Pos[pid][j]*LSCALE-p->Cen[j];
		p->vel[i][j]=Pdat.Vel[pid][j]-p->VCen[j];
	  }
	  if(Tree.Branch[BranchID].Node[0].SnapNum==FinalSnap)
		p->mass[i]=Tree.Branch[BranchID].Node[0].nP;
	  else
		p->mass[i]=0;

	  HBTInt imax=get_Mmax_Node(BranchID);
	  p->massTVV[i][0]=Tree.Branch[BranchID].Node[imax].nP;
	  p->snapTVV[i][0]=Tree.Branch[BranchID].Node[imax].SnapNum;
	  p->IsDirectInfall[i]=Tree.Branch[BranchID].Node[0].nP_infall; //overloaded with DTree infall mass
	  p->massTVV[i][1]=Tree.Branch[BranchID].Node[0].nP_infall;
	  HBTInt InfallNodeID[2];
	  get_infall_nodes(BranchID, InfallNodeID); //not according to Rvir cross, but according to subhalo rank. return the snapshot just before it becomes a satellite if possible.
	  for(j=0;j<2;j++)
	  {
		if(InfallNodeID[j]<0)
		{
		  p->snapTVV[i][j+1]=-1;
		  p->massTVV[i][j+1]=0;
		}
		else
		{
		  p->snapTVV[i][j+1]=Tree.Branch[BranchID].Node[InfallNodeID[j]].SnapNum;
		  p->massTVV[i][j+1]=Tree.Branch[BranchID].Node[InfallNodeID[j]].nP;
		}
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
	status = H5LTmake_dataset(file_id,"/DirectInfall",1,&dims[0],H5T_NATIVE_INT,p->IsDirectInfall);//overloaded with DTree infall mass
    status = H5LTmake_dataset(file_id,"/mass",1,&dims[0],H5T_NATIVE_INT,p->mass);//in particles
    status = H5LTmake_dataset(file_id,"/snapTVV",2,dims,H5T_NATIVE_INT,p->snapTVV);//SnapTidal, SnapFirstRvInfall, SnapLastRvInfall
	status = H5LTmake_dataset(file_id,"/massTVV",2,dims,H5T_NATIVE_INT,p->massTVV);//in particles
		
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[0]: (%g,%g,%g)\n",p->pos[0][0],p->pos[0][1],p->pos[0][2]);
}

static int comp_IDatInt(const void *a, const void *b)//used to sort PID in ascending order; 
{
  if((*(IDatInt *)a) > (*(IDatInt *)b))
    return +1;

  if((*(IDatInt *)a) < (*(IDatInt *)b))
    return -1;

  return 0;
}

void load_tree(struct Tree *Tree, char *infile)
{
	hid_t    file_id,group_id;
    herr_t      status;
    size_t     i, j,nread=0;
    hsize_t dims[2];
	
	size_t nload;
    int *p, *BranchOffset;
	long int *plong;
    FloatMat A;
    GenericMat B;
    char grpname[32];
    
	if(!try_readfile(infile))
	{
	  fprintf(stderr, "Error: fail to open file %s\n", infile);
	  exit(1);
	}
    printf("loading %s...\n",infile);
	
    sprintf(B.name,"/treeIndex/firstNode");
	nload=load_hdfmatrix(infile,&B,1,H5T_NATIVE_INT);
	BranchOffset=B.x;
	Tree->NBranch=B.size[0];
	Tree->Branch=malloc(sizeof(struct TreeBranch)*Tree->NBranch);
	free(B.size);
	
	sprintf(B.name,"/treeIndex/numberOfNodes");
	nload=load_hdfmatrix(infile,&B,1,H5T_NATIVE_INT);
	p=B.x;
	HBTInt totNNodes=0;
	for(i=0;i<Tree->NBranch;i++)
	{
	  Tree->Branch[i].NNode=p[i];
	  totNNodes+=p[i];
	}
	FREE_HDFMATRIX(B)
	
	Tree->Branch[0].Node=malloc(sizeof(struct TreeNode)*totNNodes);
	for(i=0;i<Tree->NBranch;i++)
	  Tree->Branch[i].Node=Tree->Branch[0].Node+BranchOffset[i];
	free(BranchOffset);
	
	sprintf(B.name, "/haloTrees/snapshotNumber");
    load_hdfmatrix(infile,&B,1,H5T_NATIVE_INT);
    if(B.size[0]!=totNNodes)
    {
      printf("Error, number of nodes mismatch %zd,%ld\n",B.size[0],(long)totNNodes);
      exit(1);
    }
    p=B.x;
	for(i=0;i<totNNodes;i++)
	  Tree->Branch[0].Node[i].SnapNum=p[i];
    FREE_HDFMATRIX(B)
	
	sprintf(B.name, "/haloTrees/isFoFCentre");
    load_hdfmatrix(infile,&B,1,H5T_NATIVE_INT);
    p=B.x;
	for(i=0;i<totNNodes;i++)
	  Tree->Branch[0].Node[i].IsFoFCentre=p[i];
    FREE_HDFMATRIX(B)
	
	sprintf(B.name, "/haloTrees/mostBoundID");
    load_hdfmatrix(infile,&B,1,H5T_NATIVE_LONG);
    plong=B.x;
	#ifdef HBTPID_RANKSTYLE
	IDatInt *PIDs=load_PIDs_Sorted();
	for(i=0;i<totNNodes;i++)
	{
	  IDatInt tmpID=plong[i];
	  Tree->Branch[0].Node[i].MBDID=(IDatInt *)bsearch(&tmpID,PIDs,NP_DM,sizeof(IDatInt),comp_IDatInt)-PIDs;
	}
	myfree(PIDs);
	#else
	for(i=0;i<totNNodes;i++)
	  Tree->Branch[0].Node[i].MBDID=plong[i];
	#endif
    FREE_HDFMATRIX(B)
		
	sprintf(B.name, "/haloTrees/particleNumber");
    load_hdfmatrix(infile,&B,1,H5T_NATIVE_INT);
    p=B.x;
	for(i=0;i<totNNodes;i++)
	  Tree->Branch[0].Node[i].nP=p[i];
    FREE_HDFMATRIX(B)
	
	sprintf(B.name, "/haloTrees/originalParticleNumber");
    load_hdfmatrix(infile,&B,1,H5T_NATIVE_INT);
    p=B.x;
	for(i=0;i<totNNodes;i++)
	  Tree->Branch[0].Node[i].nP_infall=p[i];
    FREE_HDFMATRIX(B)
	
	sprintf(A.name, "/haloTrees/position");
    load_hdfmatrixF(infile,&A,1);
	for(i=0;i<totNNodes;i++)
	  for(j=0;j<3;j++)
		Tree->Branch[0].Node[i].Pos[j]=A.x[i*3+j];
    FREE_HDFMATRIX(A)
}

void free_tree(struct Tree *Tree)
{
  if(Tree->NBranch)
	free(Tree->Branch[0].Node);
  free(Tree->Branch);
}

int get_mostmassive_branch(struct Tree *Tree)
{
  HBTInt m, i, imax;
  for(m=0,imax=0,i=0;i<Tree->NBranch;i++)
  {
	if(m<Tree->Branch[i].Node[0].nP)
	{
	  m=Tree->Branch[i].Node[0].nP;
	  imax=i;
	}
  }
  printf("Most massive branch found at snapshot %d, with %ld particles\n", Tree->Branch[imax].Node[0].SnapNum, (long)m);
  return imax;
}

HBTInt get_Mmax_Node(HBTInt BranchID)
{
  HBTInt i,imax=0,Mmax=0,m;
  struct TreeNode *Node=Tree.Branch[BranchID].Node;
  for(i=0;i<Tree.Branch[BranchID].NNode;i++)
  {
	m=Node[i].nP;
	if(m>Mmax)
	{
	  imax=i;
	  Mmax=m;
	}
  }
  return imax;
}

void get_infall_nodes(HBTInt BranchID, HBTInt InfallNodeID[2])
{
  HBTInt i;
  struct TreeNode *Node=Tree.Branch[BranchID].Node;
  InfallNodeID[0]=-1;
  InfallNodeID[1]=-1;
  for(i=Tree.Branch[BranchID].NNode-1;i>=0;i--)
  {
	  if(Node[i].IsFoFCentre) //isolated
		InfallNodeID[1]=i; //record for last-infall
	  else //satellite
	  {
		if(InfallNodeID[0]<0)//no recorded infall yet
		{
		  if(i+1<Tree.Branch[BranchID].NNode)
			InfallNodeID[0]=i+1; //previous infall
		  else
			InfallNodeID[0]=i; //current node, because previous node is invalid
		}
	  }
  }
}