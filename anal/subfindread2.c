//SUBFIND v2 (GADGET2) FORMAT
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
  
static int Ngroups,Nids,TotNgroups,NFiles;
static int *groupIDs, *SubMostBoundID;
static int Nsubhalos, *SubParentHalo, *SubOffset, *SubLen, *FirstSubOfHalo, *NsubPerHalo;
static float *SubPos, *SubVel, *SubVelDisp, *SubVmax, *SubSpin, *SubHalfMass;
static float *Halo_M_Mean200, *Halo_R_Mean200, *Halo_M_Crit200, *Halo_R_Crit200, *Halo_M_TopHat200,
*Halo_R_TopHat200;


void subread(int SnapshotNum,char *subdir)
{
  FILE *fd;
  char buf[1024];

	#ifdef SNAPLIST
	sprintf(buf, "%s/sub_tab_%04d",subdir, snaplist[SnapshotNum]);
	#else
  	sprintf(buf, "%s/sub_tab_%03d",subdir, SnapshotNum);
  	#endif
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file %s\n",buf);
      exit(1);
    }
    
  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NFiles, sizeof(int), 1, fd);
  fread(&Nsubhalos, sizeof(int), 1, fd);
    
  NsubPerHalo = malloc(Ngroups * sizeof(int));
  FirstSubOfHalo =malloc(Ngroups * sizeof(int));

  Halo_M_Mean200 = malloc(Ngroups * sizeof(float));
  Halo_R_Mean200 = malloc(Ngroups * sizeof(float));
  Halo_M_Crit200 = malloc(Ngroups * sizeof(float));
  Halo_R_Crit200 = malloc(Ngroups * sizeof(float));
  Halo_M_TopHat200 = malloc(Ngroups * sizeof(float));
  Halo_R_TopHat200 = malloc(Ngroups * sizeof(float));

  SubLen = malloc(Nsubhalos * sizeof(int));
  SubOffset = malloc(Nsubhalos * sizeof(int));
  SubParentHalo = malloc(Nsubhalos * sizeof(int));

  SubPos = malloc(Nsubhalos * 3 * sizeof(float));//center of minimum potential
  SubVel = malloc(Nsubhalos * 3 * sizeof(float));
  SubVelDisp = malloc(Nsubhalos * sizeof(float));
  SubVmax = malloc(Nsubhalos * sizeof(float));
  SubSpin = malloc(Nsubhalos * 3 * sizeof(float));
  SubMostBoundID = malloc(Nsubhalos * sizeof(int));
  SubHalfMass = malloc(Nsubhalos * sizeof(float));

  fread(NsubPerHalo, sizeof(int), Ngroups, fd);
  fread(FirstSubOfHalo, sizeof(int), Ngroups, fd);

  fread(SubLen, sizeof(int), Nsubhalos, fd);
  fread(SubOffset, sizeof(int), Nsubhalos, fd);
  fread(SubParentHalo, sizeof(int), Nsubhalos, fd);

  fread(Halo_M_Mean200, sizeof(float), Ngroups, fd);
  fread(Halo_R_Mean200, sizeof(float), Ngroups, fd);
  fread(Halo_M_Crit200, sizeof(float), Ngroups, fd);
  fread(Halo_R_Crit200, sizeof(float), Ngroups, fd);
  fread(Halo_M_TopHat200, sizeof(float), Ngroups, fd);
  fread(Halo_R_TopHat200, sizeof(float), Ngroups, fd);


  fread(SubPos, 3 * sizeof(float), Nsubhalos, fd);
  fread(SubVel, 3 * sizeof(float), Nsubhalos, fd);
  fread(SubVelDisp, sizeof(float), Nsubhalos, fd);
  fread(SubVmax, sizeof(float), Nsubhalos, fd);
  fread(SubSpin, 3 * sizeof(float), Nsubhalos, fd);
  fread(SubMostBoundID, sizeof(int), Nsubhalos, fd);
  fread(SubHalfMass, sizeof(float), Nsubhalos, fd);

  fclose(fd);
  printf("tab ok\n");fflush(stdout);
  
  #ifdef SNAPLIST
  sprintf(buf, "%s/sub_ids_%04d",subdir, snaplist[SnapshotNum]);
  #else
  sprintf(buf, "%s/sub_ids_%03d",subdir, SnapshotNum);
  #endif
  if(!(fd = fopen(buf, "r")))
    {
      printf( "can't open file \n");
      exit(1);
    }
    
  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  printf("Nids:%d\n",Nids);fflush(stdout);  
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NFiles, sizeof(int), 1, fd);

  groupIDs = malloc(sizeof(int) * Nids);
  
  fread(groupIDs, sizeof(int), Nids, fd);

  fclose(fd);
  }


#ifdef PARAM_FILE_INCLUDED
void subfindcat2BT(SUBCATALOGUE *SubCat)
{
	int subid,grpid,i,j;
	
	SubCat->Ngroups=Ngroups;
	SubCat->Nsubs=Nsubhalos;
	SubCat->Nids=SubOffset[Nsubhalos-1]+SubLen[Nsubhalos-1];
	SubCat->GrpLen_Sub=NsubPerHalo;
	SubCat->GrpOffset_Sub=FirstSubOfHalo;
	if(FirstSubOfHalo[0]!=0)
	{
	printf("error,subhaloid not start from 0\n");
	exit(1);
	}
	if(SubParentHalo[0]!=0)
	{
	printf("error,Host haloid not start from 0\n");
	exit(1);
	}
	SubCat->SubLen=SubLen;
	SubCat->SubOffset=SubOffset;
	
	SubCat->HaloChains=mymalloc(sizeof(struct Chain_data)*Nsubhalos);
	SubCat->SubRank=mymalloc(sizeof(int)*Nsubhalos);
	SubCat->Property=mymalloc(sizeof(struct SubProperty)*Nsubhalos);
	SubCat->PSubArr=mymalloc(sizeof(int *)*Nsubhalos);
	for(grpid=0;grpid<Ngroups;grpid++)
	{
		for(i=0,subid=FirstSubOfHalo[grpid];subid<FirstSubOfHalo[grpid]+NsubPerHalo[grpid];subid++,i++)
		{
			SubCat->SubRank[subid]=i;
			SubCat->PSubArr[subid]=groupIDs+SubOffset[subid];
			SubCat->HaloChains[subid].HostID=SubParentHalo[subid];
			for(j=0;j<3;j++)
			{
			SubCat->Property[subid].CoM[j]=SubPos[subid*3+j];//center of minimum potential
			SubCat->Property[subid].VCoM[j]=SubVel[subid*3+j];//average vel?
			}
		}
	}
	//~ for(subid=0;subid<Nsubhalos-1;subid++)
	//~ {
		//~ if(SubLen[subid]+SubOffset[subid]!=SubOffset[subid+1])
		//~ {
		//~ printf("error,lastsubid=%d,suboffset=%d,sub.offset=%d,lastsublen=%d,lastsuboffset=%d\n",subid,SubOffset[subid+1],SubLen[subid]+SubOffset[subid],SubLen[subid],SubOffset[subid]);
		//~ exit(1);
		//~ }
	//~ }
	SubCat->sub_hierarchy=NULL;
	SubCat->Nbirth=0;
	SubCat->Ndeath=0;
	SubCat->NQuasi=0;
	SubCat->Nsplitter=0;
}
void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir)
{
	subread(Nsnap,inputdir);
	subfindcat2BT(SubCat);
}
#ifdef STANDALONE
int main()
{
	int Nsnap=86;
	char inputdir[1024]="/gpfs/data/arj/projects/jiaxin/run";
	subread(Nsnap,inputdir);
	return 0;
}
#endif
#endif	
