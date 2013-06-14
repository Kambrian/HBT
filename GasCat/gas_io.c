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

void load_gas_data(HBTInt Nsnap, char *SnapPath)
{
HBTInt i,j,Ngas,Ndm,Nother,Nmass,dummy,dummy2;
	#ifdef SKIP
		#undef SKIP
		#undef SKIP2
		#undef CHECK
	#endif
	#define SKIP fread(&dummy,sizeof(dummy),1,fp)
	#define SKIP2 fread(&dummy2,sizeof(dummy2),1,fp)
	#define CHECK if(dummy!=dummy2){fprintf(logfile,"error!record brackets not match for Snap%d!\t%d,%d\n",Nsnap,dummy,dummy2);fflush(logfile);exit(1);} 
	FILE *fp;
	char buf[1024];

struct
{
IDatInt *PID;
IDatReal (*Pos)[3];
IDatReal (*Vel)[3];
IDatReal *U;
IDatReal *Rho;
} IDat;

	sprintf(buf,"%s/%s_%03d",SnapPath,SNAPFILE_BASE,Nsnap);
	if((fp=fopen(buf,"r"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		fflush(logfile);
		exit(1);
	}
	SKIP;
	fread(&header,sizeof(header),1,fp);
	SKIP2;
	CHECK;

	header.Hz=HUBBLE0 * sqrt(header.Omega0 / (header.time * header.time * header.time) 
			+ (1 - header.Omega0 - header.OmegaLambda) / (header.time * header.time)
			+ header.OmegaLambda);//Hubble param for the current catalogue;
	
	Ngas=header.npart[0];	Ndm=header.npart[1];	Nother=header.npart[2]+header.npart[3]+header.npart[4]+header.npart[5];
	if(Ngas<=0)
	{
		fprintf(stderr,"warning,no gas particles:%d\n",Ngas);
		return;
	}
	for(Nmass=0,i=0;i<6;i++)
	{
		if(!header.mass[i]) Nmass+=header.npart[i];
	}

	IDat.Pos=mymalloc(sizeof(IDatReal)*3*Ngas);
	SKIP;
	fread(IDat.Pos,sizeof(IDatReal)*3,Ngas,fp);
	fseek(fp,(Ndm+Nother)*sizeof(IDatReal)*3,SEEK_CUR);
	SKIP2;
	CHECK;
	
	IDat.Vel=mymalloc(sizeof(IDatReal)*3*Ngas);
	SKIP;
	fread(IDat.Vel,sizeof(IDatReal)*3,Ngas,fp);
	fseek(fp,(Ndm+Nother)*sizeof(IDatReal)*3,SEEK_CUR);
	SKIP2;
	CHECK;
	
	IDat.PID=mymalloc(sizeof(IDatInt)*Ngas);
	SKIP;
	fread(IDat.PID,sizeof(IDatInt),Ngas,fp);
	fseek(fp,(Ndm+Nother)*sizeof(IDatInt),SEEK_CUR);
	SKIP2;
	CHECK;
	
	if(Nmass)
	{
		SKIP;
		fseek(fp,Nmass*sizeof(IDatReal),SEEK_CUR);
		SKIP2;
		CHECK;
	}
	
	IDat.U=mymalloc(sizeof(IDatReal)*Ngas);
	SKIP;
	fread(IDat.U,sizeof(IDatReal),Ngas,fp);
	SKIP2;
	CHECK;

	IDat.Rho=mymalloc(sizeof(IDatReal)*Ngas);
	SKIP;
	fread(IDat.Rho,sizeof(IDatReal),Ngas,fp);
	SKIP2;
	CHECK;
	
	fclose(fp);
	#undef SKIP
	#undef SKIP2
	#undef CHECK

#ifdef HBTPID_RANKSTYLE
  #error HBTPID_RANKSTYLE not yet implemented for GasData
#endif
	
#ifdef SAME_INTTYPE
  Gdat.PID=IDat.PID;
#else
  Gdat.PID=mymalloc(sizeof(HBTInt)*Ngas);
  for(i=0;i<Ngas;i++)
	Gdat.PID[i]=IDat.PID[i];
  myfree(IDat.PID);
#endif

#ifdef SAME_REALTYPE
	  Gdat.Pos=IDat.Pos;
#else
	  Gdat.Pos=mymalloc(sizeof(HBTxyz)*NP_GAS);
	  for(i=0;i<NP_GAS;i++)
		for(j=0;j<3;j++)
			Gdat.Pos[i][j]=IDat.Pos[i][j];
	  myfree(IDat.Pos);
#endif
  	
  #ifdef CONVERT_LENGTH_MPC_KPC
  for(i=0;i<NP_GAS;i++)
  for(j=0;j<3;j++)
  Gdat.Pos[i][j]*=1000.;
  #endif
  
  #ifdef PERIODIC_BDR
  for(i=0;i<NP_GAS;i++)
  for(j=0;j<3;j++)
  Gdat.Pos[i][j]=position_modulus(Gdat.Pos[i][j]);
  #endif	
  
#ifdef SAME_REALTYPE
	  Gdat.Vel=IDat.Vel;
#else
	  Gdat.Vel=mymalloc(sizeof(HBTxyz)*NP_GAS);
	  for(i=0;i<NP_GAS;i++)
		for(j=0;j<3;j++)
			Gdat.Vel[i][j]=IDat.Vel[i][j];
	  myfree(IDat.Vel);
#endif

#ifdef SAME_REALTYPE
	  Gdat.U=IDat.U;
	  Gdat.Rho=IDat.Rho;
#else
	  Gdat.U=mymalloc(sizeof(HBTReal)*NP_GAS);
	  Gdat.Rho=mymalloc(sizeof(HBTReal)*NP_GAS);
	  for(i=0;i<NP_GAS;i++)
	  {
		Gdat.U[i]=IDat.U[i];
		Gdat.Rho[i]=IDat.Rho[i];
	  }
	  myfree(IDat.U);
	  myfree(IDat.Rho);
#endif
	
}
void free_gas_data()
{
	myfree(Gdat.PID);
	myfree(Gdat.Pos);
	myfree(Gdat.Vel);
	myfree(Gdat.U);
	myfree(Gdat.Rho);
}
void save_gashalocat(HBTInt Nsnap, GASHALOCAT *Cat,char *gaspath)
{
FILE *fd;
char buf[1024];
HBTInt NFiles;
NFiles=1;

  sprintf(buf, "%s/gas_group_%03d",gaspath,Nsnap);
  if(!(fd = fopen(buf, "w")))
    {
	printf("can't open file `%s'\n", buf);
	exit(1);
    }

  fwrite(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nids, sizeof(HBTInt), 1, fd);

  fwrite(Cat->Len, sizeof(HBTInt), Cat->Ngroups, fd);
  fwrite(Cat->Offset,sizeof(HBTInt), Cat->Ngroups, fd);

  fwrite(Cat->PIDorIndex, sizeof(HBTInt), Cat->Nids, fd);
  fclose(fd);
}
void load_gashalocat(HBTInt Nsnap,GASHALOCAT *Cat,char *GrpPath)
{
  FILE *fd;
  HBTInt i;
  char buf[1024];

  sprintf(buf, "%s/gas_group_%03d",GrpPath,Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);

  fprintf(logfile,"Snap=%d  Ngroups=%d  Nids=%d\n", Nsnap, Cat->Ngroups, Cat->Nids);fflush(logfile);

  Cat->Len= mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
  //~ Cat->HaloMask=mymalloc(sizeof(short)*NP_DM);
  //~ Cat->HaloMaskSrc=mymalloc(sizeof(short)*NP_DM);
  //~ Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_GAS);//.............this can be discarded;
  Cat->ID2Halo=NULL;//different from load_group_cat

  fread(Cat->Len, sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->Offset,sizeof(HBTInt), Cat->Ngroups, fd);

  fread(Cat->PIDorIndex, sizeof(HBTInt), Cat->Nids, fd);
	for(i=0;i<Cat->Nids;i++)
	{
 	  if(Cat->PIDorIndex[i] == 0)
		fprintf(logfile,"i=%d gasPID=%d\n", i, (HBTInt)Cat->PIDorIndex[i]);//check if PID begin with ID=1;
  	  }

  fclose(fd);
}

void free_gashalocat(GASHALOCAT *A)
{
	free(A->Len);
	free(A->Offset);
	free(A->PIDorIndex);
	//~ free(A->HaloMask);
	//~ free(A->HaloMaskSrc);
	//~ free(A->ID2Halo);//different from free_group_cat
}

HBTInt prepare_gasind2halo(GASHALOCAT *A)
{HBTInt i,haloid,pid;
//	A->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);

	#pragma omp parallel 
		{
	#pragma omp for
	for(i=0;i<NP_GAS;i++)//initialization
	{
		A->ID2Halo[i]=-1;/*return -1 if the PID does not belong to a Halo,
								i.e,we consider the backgroud as a halo with haloid=-1; 
								note that this only make sense when we try to find host for bound structures */
	}
	#pragma omp for private(i,pid,haloid)
	for(haloid=0;haloid<A->Ngroups;haloid++)
	{
		for(i=0;i<A->Len[haloid];i++)
		{
			pid=A->PIDorIndex[A->Offset[haloid]+i];//Pindex ranges [0,NP_DM);
			A->ID2Halo[pid]=haloid;//haloIDs begins from id=0
		}
	}
		}
	return 1; // means id2halo ready;
}

void fresh_gasID2Index(const void *src, HBTInt src_len)
{
	HBTInt i,j,pid,*Arr;			//src_len==ArrLen;
	GASHALOCAT *Cat;			//src_len==-1,or FRSH_GRPCAT
	GASSUBCAT *SubCat;			//src_len==-2,or FRSH_SUBCAT
	GASSRCCAT *SrcCat;			//src_len==-3,or FRSH_SRCCAT
	
	#ifndef PID_ORDERED
	HBTInt *PIndex;
	PIndex=calloc(NP_SIM,sizeof(HBTInt));//alloc and fill with NULL
	PIndex=PIndex-1;/*shift PIndex by -1 element so that it is accessed through PIndex[1~N],
											i.e.,its subscript ranges the same as particle ID range.*/	
									//since PID begins from ID=1, PIndex[0] will not be used;PID_all[PIndex[PID]]=PID;
	/*====make ID index for query====*/
	for(i=0;i<NP_GAS;i++)
		PIndex[Gdat.PID[i]]=i;		//suppose PID begins from ID=1;
	#endif
	
	if(FRSH_GRPCAT==src_len)//cat
	{
		#ifndef GRPINPUT_INDEX
		Cat=(GASHALOCAT *)src;
		for(i=0;i<Cat->Nids;i++)
			#ifdef PID_ORDERED
			Cat->PIDorIndex[i]--;//change from ID [1,N] to Index [0,N-1]
			#else
			Cat->PIDorIndex[i]=PIndex[Cat->PIDorIndex[i]];//Now PIDorIndex has been refilled with particle Index in Gdat
			#endif
		#endif	
	}
	else if(FRSH_SUBCAT==src_len)//subcat
	{
		SubCat=(GASSUBCAT *)src;
		/*====refresh SubCat with particle Index===*/
		for(i=0;i<SubCat->Nsubs;i++)
		{
			for(j=0;j<SubCat->SubLen[i];j++)
			{
				pid=SubCat->PSubArr[i][j];
				if(pid<1||pid>NP_SIM)
				{
					fprintf(logfile,"(%d,%d),%d,%d\n",i,j,pid,SubCat->SubLen[i]);fflush(logfile);
					exit(1);
				}
				#ifdef PID_ORDERED
				SubCat->PSubArr[i][j]--;
				#else
				SubCat->PSubArr[i][j]=PIndex[pid];
				#endif
			}
		}
	}
	else if(FRSH_SRCCAT==src_len)//srccat
	{
		SrcCat=(GASSRCCAT *)src;
		/*====refresh SubCat with particle Index===*/
		for(i=0;i<SrcCat->Nsubs;i++)
		{
			for(j=0;j<SrcCat->SubLen[i];j++)
			{
				pid=SrcCat->PSubArr[i][j];
				if(pid<1||pid>NP_SIM)
				{
					fprintf(logfile,"(%d,%d),%d,%d\n",i,j,pid,SrcCat->SubLen[i]);fflush(logfile);
					exit(1);
				}
				#ifdef PID_ORDERED
				SrcCat->PSubArr[i][j]--;
				#else
				SrcCat->PSubArr[i][j]=PIndex[pid];
				#endif
			}
		}
	}
	else if(src_len>=0)
	{
		Arr=(HBTInt *)src;
		for(i=0;i<src_len;i++)
		#ifdef PID_ORDERED
			Arr[i]--;
		#else
			Arr[i]=PIndex[Arr[i]];
		#endif
	}
	else
	{
		 fprintf(logfile,"Error: wrong srclen when using fresh_gasID2Index\nsrclen must be an integer no smaller than -3\n");fflush(logfile);
		exit(1);
	}
	#ifndef PID_ORDERED
		free(PIndex+1);
	#endif
}

void save_gassrccat(HBTInt Nsnap,GASSRCCAT *Cat,char *outputdir)
{
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/gassrccat_%03d", outputdir, Nsnap);
  myfopen(fd,buf,"w");

  fwrite(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fwrite(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fwrite(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
for(i=0;i<Cat->Nsubs;i++)
  fwrite(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);

  fclose(fd);
}	
void load_gassrccat(HBTInt Nsnap,GASSRCCAT *Cat,char *inputdir)
{
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/gassrccat_%03d", inputdir, Nsnap);
  myfopen(fd,buf,"r");

  fread(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  create_gassrccat(Cat);
  fread(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
for(i=0;i<Cat->Nsubs;i++)
{
  Cat->PSubArr[i]=mymalloc(sizeof(HBTInt)*Cat->SubLen[i]);
  fread(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
}

  fclose(fd);
}
void save_gassubcat(HBTInt Nsnap,GASSUBCAT *Cat,char *outputdir)
{
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/gassubcat_%03d", outputdir, Nsnap);
  myfopen(fd,buf,"w");

  fwrite(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fwrite(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fwrite(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
for(i=0;i<Cat->Nsubs;i++)
  fwrite(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);

  fwrite(Cat->Property,sizeof(struct GasProperty),Cat->Nsubs,fd);
  fclose(fd);
}
void load_gassubcat(HBTInt Nsnap,GASSUBCAT *Cat,char *outputdir)
{
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/gassubcat_%03d", outputdir, Nsnap);
  myfopen(fd,buf,"r");

  fread(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  create_gassubcat(Cat);
  fread(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
for(i=0;i<Cat->Nsubs;i++)
{
  Cat->PSubArr[i]=mymalloc(sizeof(HBTInt)*Cat->SubLen[i]);
  fread(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
}

  fread(Cat->Property,sizeof(struct GasProperty),Cat->Nsubs,fd);
  fclose(fd);
}
void create_gassrccat(GASSRCCAT *GSrcCat)
{
	HBTInt i,Nsubs;
	Nsubs=GSrcCat->Nsubs;//this must be assigned before using create()
	GSrcCat->SubLen=mymalloc(sizeof(HBTInt)*Nsubs);
	GSrcCat->SubOffset=mymalloc(sizeof(HBTInt)*Nsubs);
	GSrcCat->PSubArr=mymalloc(sizeof(HBTInt *)*Nsubs);
	for(i=0;i<Nsubs;i++)
		GSrcCat->PSubArr[i]=NULL;
}
void free_gassrccat(GASSRCCAT *GSrcCat)
{
	myfree(GSrcCat->SubLen);
	myfree(GSrcCat->SubOffset);
	myfree(GSrcCat->PSubArr);
}
void erase_gassrccat(GASSRCCAT *GSrcCat)
{
	HBTInt i;
	for(i=0;i<GSrcCat->Nsubs;i++)
		myfree(GSrcCat->PSubArr[i]);
	free_gassrccat(GSrcCat);
}
void create_gassubcat(GASSUBCAT *GSubCat)
{
	HBTInt i,Nsubs;
	Nsubs=GSubCat->Nsubs;//this must be assigned before using create()
	GSubCat->SubLen=mymalloc(sizeof(HBTInt)*Nsubs);
	GSubCat->SubOffset=mymalloc(sizeof(HBTInt)*Nsubs);
	GSubCat->PSubArr=mymalloc(sizeof(HBTInt *)*Nsubs);
	for(i=0;i<Nsubs;i++)
		GSubCat->PSubArr[i]=NULL;
	GSubCat->Property=mymalloc(sizeof(struct GasProperty)*Nsubs);
}
void free_gassubcat(GASSUBCAT *GSubCat)
{
	myfree(GSubCat->SubLen);
	myfree(GSubCat->SubOffset);
	myfree(GSubCat->PSubArr);
	myfree(GSubCat->Property);
}
void erase_gassubcat(GASSUBCAT *GSubCat)
{
	HBTInt i;
	for(i=0;i<GSubCat->Nsubs;i++)
		myfree(GSubCat->PSubArr[i]);
	free_gassubcat(GSubCat);
}
void complete_N_save_gas(GASSRCCAT *GSrcCat,GASSUBCAT *GSubCat,HBTInt SnapshotNum,char *gasdir)
{
	HBTInt i,j,Offset;
	
	Offset=0;
	for(i=0;i<GSrcCat->Nsubs;i++)
	{
		GSrcCat->SubOffset[i]=Offset;
		Offset+=GSrcCat->SubLen[i];
		for(j=0;j<GSrcCat->SubLen[i];j++)
			GSrcCat->PSubArr[i][j]=Gdat.PID[GSrcCat->PSubArr[i][j]];
	}
	GSrcCat->Nids=Offset;
	save_gassrccat(SnapshotNum,GSrcCat,gasdir);
	
	Offset=0;
	for(i=0;i<GSubCat->Nsubs;i++)
	{
		GSubCat->SubOffset[i]=Offset;
		Offset+=GSubCat->SubLen[i];
		for(j=0;j<GSubCat->SubLen[i];j++)
			GSubCat->PSubArr[i][j]=Gdat.PID[GSubCat->PSubArr[i][j]];
	}
	GSubCat->Nids=Offset;
	save_gassubcat(SnapshotNum,GSubCat,gasdir);
}
