/*to add SubProperty to old subcatalogues (version <= 7.4);
 * saved in v7.6 format (fixed the bug of not restoring Ind2ID before saving)
 * Difference between results applying this patch and those from BT76:
 * 		VCoM from CoreFrac0 percent Most-Bound particles       for patched version
 * 			 from CoreFrac  percent Lowest-Potential particles for BT76*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

void load_sub_catalogue_74(int Nsnap, SUBCATALOGUE *Cat, char *SubCatPath);

int main(int argc,char **argv)
{
char outputdir[1024];
SUBCATALOGUE SubCat;
int i,j,pid,SnapshotNum,subid,CoreLen,SubLen,*SubPIndex;
float *CoM;
double pot,SubPot,SubKin;
double Vcomx,Vcomy,Vcomz,SubAMx,SubAMy,SubAMz;
double Hz,sqa,Time,PartMass;
double dx,dy,dz,dvx,dvy,dvz;
int snapbegin,snapend;
char buf[1024];

if(argc!=3)
{
	printf("usage: %s [snap_begin] [snap_end]\n",argv[0]);
	return 1;
}
snapbegin=atoi(argv[1]);
snapend=atoi(argv[2]);

sprintf(outputdir,"%s/patched",SUBCAT_DIR);
sprintf(buf,"%s/logfile.patch",outputdir);
myfopen(logfile,buf,"w");
for(SnapshotNum=snapbegin;SnapshotNum<=snapend;SnapshotNum++)
{
	load_sub_catalogue_74(SnapshotNum,&SubCat,SUBCAT_DIR);
	load_particle_data(SnapshotNum,SNAPSHOT_DIR);
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	Hz=header.Hz;
	PartMass=header.mass[1];	
	fresh_ID2Index(&SubCat,FRSH_SUBCAT);
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		SubLen=SubCat.SubLen[subid];
		if(0==SubLen)
		{
		SubCat.Property[subid].VCoM[0]=SubCat.Property[subid].VCoM[1]=SubCat.Property[subid].VCoM[2]=0.;
		SubCat.Property[subid].Pot=0.;
		SubCat.Property[subid].Kin=0.;
		SubCat.Property[subid].AM[0]=SubCat.Property[subid].AM[0]=SubCat.Property[subid].AM[0]=0.;			
		}
		else
		{
		SubPIndex=SubCat.PSubArr[subid];
		CoreLen=SubLen*CoreFrac0;
		CoreLen=(CoreLen<CoreLenMin)?CoreLenMin:CoreLen;
		Vcomx=Vcomy=Vcomz=0.;
		#pragma omp parallel for private(i,pid) reduction(+:Vcomx,Vcomy,Vcomz) \
								 schedule(dynamic,1) if(CoreLen>NParaMin) 
		for(i=0;i<CoreLen;i++)
		{
			pid=SubPIndex[i];
			Vcomx+=Pdat.Vel[pid][0];
			Vcomy+=Pdat.Vel[pid][1];
			Vcomz+=Pdat.Vel[pid][2];
		}
		Vcomx/=CoreLen;Vcomy/=CoreLen;Vcomz/=CoreLen;
		SubCat.Property[subid].VCoM[0]=Vcomx*sqa;//physical vel
		SubCat.Property[subid].VCoM[1]=Vcomy*sqa;
		SubCat.Property[subid].VCoM[2]=Vcomz*sqa;
		
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubPIndex,Pdat.Pos);
		CoM=SubCat.Property[subid].CoM;
		SubPot=SubKin=SubAMx=SubAMy=SubAMz=0.;
		#pragma omp parallel for private(i,pot,dvx,dvy,dvz,dx,dy,dz) reduction(+:SubPot,SubKin,SubAMx,SubAMy,SubAMz) \
								schedule(dynamic,1) if(SubLen>NParaMin) 
		for(i=0;i<SubLen;i++)
		{
			pot=tree_treeevaluate_potential(Pdat.Pos[SubPIndex[i]],SubPIndex,Pdat.Pos);
			pot=(PartMass*pot+PartMass/SofteningHalo)*G/Time;/*exclude self-potential 
															*which was included when evaluating potential
															*  (-2.8M/h=-M/softening when r=0)*/
			SubPot+=pot;			
																
			 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-Vcomx);//relative vel.
			 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-Vcomy);
			 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-Vcomz);
			 dx=Pdat.Pos[SubPIndex[i]][0]-CoM[0];
			 dy=Pdat.Pos[SubPIndex[i]][1]-CoM[1];
			 dz=Pdat.Pos[SubPIndex[i]][2]-CoM[2];
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
		SubCat.Property[subid].Pot=SubPot/SubLen;	
		SubCat.Property[subid].Kin=SubKin/SubLen;
		SubCat.Property[subid].AM[0]=SubAMx/SubLen;
		SubCat.Property[subid].AM[1]=SubAMy/SubLen;
		SubCat.Property[subid].AM[2]=SubAMz/SubLen;
		}
	}
	for(subid=0;subid<SubCat.Nsubs;subid++)
		for(i=0;i<SubCat.SubLen[subid];i++)
			SubCat.PSubArr[subid][i]=Pdat.PID[SubCat.PSubArr[subid][i]];
	save_sub_catalogue(SnapshotNum,&SubCat,outputdir);
	erase_sub_catalogue(&SubCat);
}
fclose(logfile);
return 0;
}

void load_sub_catalogue_74(int Nsnap, SUBCATALOGUE *Cat, char *SubCatPath)
{
FILE *fd;
char buf[1024];
int i;

  sprintf(buf, "%s/subcat_%03d", SubCatPath, Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Ngroups, sizeof(int), 1, fd);
        Cat->GrpOffset_Sub=mymalloc(sizeof(int)*Cat->Ngroups);   
	Cat->GrpLen_Sub=mymalloc(sizeof(int)*Cat->Ngroups);
  fread(&Cat->Nsubs, sizeof(int), 1, fd);
	Cat->SubLen=mymalloc(sizeof(int)*Cat->Nsubs);
	Cat->SubOffset=mymalloc(sizeof(int)*Cat->Nsubs);
	Cat->SubRank=mymalloc(sizeof(int)*Cat->Nsubs);
	Cat->HaloChains=mymalloc(sizeof(struct Chain_data)*Cat->Nsubs);
	Cat->Property=mymalloc(sizeof(struct SubProperty)*Cat->Nsubs);
	Cat->sub_hierarchy=mymalloc(sizeof(struct Hierarchy)*Cat->Nsubs);
	Cat->PSubArr=mymalloc(sizeof(int *)*Cat->Nsubs);
  fread(&Cat->Nids, sizeof(int), 1, fd);
  fread(Cat->GrpLen_Sub, sizeof(int), Cat->Ngroups, fd);
  fread(Cat->GrpOffset_Sub,sizeof(int), Cat->Ngroups, fd);
  fread(Cat->SubLen, sizeof(int), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(int), Cat->Nsubs, fd);
  fread(Cat->SubRank,sizeof(int), Cat->Nsubs, fd);
  fread(Cat->HaloChains,sizeof(struct Chain_data), Cat->Nsubs, fd);
  for(i=0;i<Cat->Nsubs;i++)
    fread(Cat->Property[i].CoM,sizeof(float), 3, fd);
  fread(Cat->sub_hierarchy,sizeof(struct Hierarchy),Cat->Nsubs,fd);
for(i=0;i<Cat->Nsubs;i++)
{
	Cat->PSubArr[i]=mymalloc(sizeof(int)*Cat->SubLen[i]);
	fread(Cat->PSubArr[i], sizeof(int), Cat->SubLen[i], fd);
}
fread(&Cat->Nbirth,sizeof(int),1,fd);
fread(&Cat->NQuasi,sizeof(int),1,fd);
fread(&Cat->Ndeath,sizeof(int),1,fd);
fread(&Cat->Nsplitter,sizeof(int),1,fd);

  fclose(fd);
}
