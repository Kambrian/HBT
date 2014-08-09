#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"


void complete_N_save(SUBCATALOGUE *SubCatB,SRCCATALOGUE *SrcCatB,HBTInt SnapshotNum,char *outputdir)
{
		HBTInt i,SubOffset,subhaloid;
	/*============fill up subcat offset======================*/			
	SubOffset=0;
	for(i=0;i<SubCatB->Nsubs;i++)
	{	
		SubCatB->SubOffset[i]=SubOffset;
		SubOffset+=SubCatB->SubLen[i];
	}	
	SubCatB->Nids=SubOffset;
		/*============fill up srccat offset======================*/			
	SubOffset=0;
	for(i=0;i<SrcCatB->Nsubs;i++)
	{	
		SrcCatB->SubOffset[i]=SubOffset;
		SubOffset+=SrcCatB->SubLen[i];
	}	
	SrcCatB->Nids=SubOffset;
	/*==============restore_subcatalogue_pid==============*/
	#ifdef PID_ORDERED
	#pragma omp parallel for private(subhaloid,i) //schedule(dynamic)
	for(subhaloid=0;subhaloid<SubCatB->Nsubs;subhaloid++)//change from index [0,NP_DM-1] to ID [1,NP_DM]
	{
		for(i=0;i<SubCatB->SubLen[subhaloid];i++)
			SubCatB->PSubArr[subhaloid][i]++;
		for(i=0;i<SrcCatB->SubLen[subhaloid];i++)
			SrcCatB->PSubArr[subhaloid][i]++;
		for(i=0;i<SrcCatB->SubLen2[subhaloid];i++)
			SrcCatB->PSubArr2[subhaloid][i]++;
	}
	#else
	#pragma omp parallel for private(subhaloid,i) //schedule(dynamic)
	for(subhaloid=0;subhaloid<SubCatB->Nsubs;subhaloid++)
	{
		for(i=0;i<SubCatB->SubLen[subhaloid];i++)
			SubCatB->PSubArr[subhaloid][i]=Pdat.PID[SubCatB->PSubArr[subhaloid][i]];
		for(i=0;i<SrcCatB->SubLen[subhaloid];i++)
			SrcCatB->PSubArr[subhaloid][i]=Pdat.PID[SrcCatB->PSubArr[subhaloid][i]];
		for(i=0;i<SrcCatB->SubLen2[subhaloid];i++)
			SrcCatB->PSubArr2[subhaloid][i]=Pdat.PID[SrcCatB->PSubArr2[subhaloid][i]];
	}
	#endif
	#pragma omp parallel sections
	{
	#pragma omp section	
	save_sub_catalogue(SnapshotNum,SubCatB,outputdir);
	#pragma omp section
	save_src_catalogue(SnapshotNum,SrcCatB,outputdir);
	}
}

//=====general ID2Index table======//
/* the hash-table implementation here is by sorting PID and use binsearch to locate keys*/
/* more general functions hcreate(),hsearch()... exists in glibc; but binsearch should be
 * more efficient than the general hsearch() I guess?*/
static int comp_PIDArr(const void *a, const void *b)//used to sort PID in ascending order; 
{
  if(((struct ID2Ind *) a)->PID > ((struct ID2Ind *) b)->PID)
    return +1;

  if(((struct ID2Ind *) a)->PID < ((struct ID2Ind *) b)->PID)
    return -1;

  return 0;
}
void fill_PIDHash2()
{
	HBTInt i;
	PIDHash.table=mymalloc(sizeof(struct ID2Ind)*NP_DM);
	for(i=0;i<NP_DM;i++)
	{
		#ifdef PID_ORDERED
		PIDHash.table[i].PID=i+1;
		#else
		PIDHash.table[i].PID=Pdat.PID[i];
		#endif
		PIDHash.table[i].PInd=i;
	}
	PIDHash.np=NP_DM;
	qsort(PIDHash.table,PIDHash.np,sizeof(struct ID2Ind),comp_PIDArr);
}
void free_PIDHash2()
{
	if(PIDHash.np)
		PIDHash.np=0;
	else
		fprintf(logfile,"warning: attempting to free PIDHash which doesn't seem to have been filled\n");
	myfree(PIDHash.table);
}
static int comp_PIDKey(const void *a, const void *b)//used to sort PID in ascending order; 
{
  if((*(HBTInt *)a) > ((struct ID2Ind *) b)->PID)
    return +1;

  if((*(HBTInt *)a) < ((struct ID2Ind *) b)->PID)
    return -1;

  return 0;
}
HBTInt lookup_PIDHash2(HBTInt PID)
{
	if(PID<0) return -1;
	
	struct ID2Ind *p;
	p=bsearch(&PID,PIDHash.table,PIDHash.np,sizeof(struct ID2Ind),comp_PIDKey);
	
	if(NULL==p)	return -1;  //no match
	
	return p->PInd;
}
//==========simple ID2Index table===========//
void fill_PIDHash1()
{
	#ifndef PID_ORDERED
	HBTInt i,idmax,idmin;
	idmax=Pdat.PID[0];idmin=Pdat.PID[0];
	for(i=1;i<NP_DM;i++)
	{
		if(Pdat.PID[i]>idmax)
			idmax=Pdat.PID[i];
		else if(Pdat.PID[i]<idmin)	
			idmin=Pdat.PID[i];
	}
	PIDHash.PIndex=mymalloc(sizeof(HBTInt)*(idmax-idmin+1));
	PIDHash.PIndex-=idmin;/*shift PIndex by idmin element so that it is accessed through PIndex[idmin~idmax],
							i.e.,its subscript ranges the same as particle ID range.
							* PID_all[PIndex[PID]]=PID;*/
	PIDHash.PInd_offset=idmin;
	for(i=idmin;i<=idmax;i++)
		PIDHash.PIndex[i]=-1; //initialize with -1, although not necessary
	/*====make ID index for query====*/
	for(i=0;i<NP_DM;i++)
	PIDHash.PIndex[Pdat.PID[i]]=i;		
	PIDHash.np=idmax-idmin+1;
	#else
	PIDHash.np=NP_DM;
	#endif
}
void free_PIDHash1()
{
	if(PIDHash.np)
		PIDHash.np=0;
	else
		fprintf(logfile,"warning: attempting to free PIDHash which doesn't seem to have been filled\n");
	#ifndef PID_ORDERED
	myfree(PIDHash.PIndex+PIDHash.PInd_offset);
	#endif
}
HBTInt lookup_PIDHash1(HBTInt PID)
{
	#ifdef PID_ORDERED
	return PID-1; //change from ID [1,N] to Index [0,N-1]
	#else
	if(PID<0)//no match
		return -1;
	return PIDHash.PIndex[PID];
	#endif	
}
//========wrapper==========//
void fill_PIDHash()
{
	#if defined PID_NEED_HASH && !defined HBTPID_RANKSTYLE
	fill_PIDHash2();
	#else
	fill_PIDHash1();
	#endif
}
void free_PIDHash()
{
	#if defined PID_NEED_HASH && !defined HBTPID_RANKSTYLE
	free_PIDHash2();
	#else
	free_PIDHash1();
	#endif
}

void fresh_ID2Index(const void *src, HBTInt src_len)
{
	HBTInt i,j,pid,*Arr;				//src_len==ArrLen;
	CATALOGUE *Cat;					//src_len==-1,or FRSH_GRPCAT
	SUBCATALOGUE *SubCat;			//src_len==-2,or FRSH_SUBCAT
	SRCCATALOGUE *SrcCat;			//src_len==-3,or FRSH_SRCCAT
//	MBDCATALOGUE ＊MbdCat;          //src_len==-4,or FRSH_MBDCAT
	
	if(0==PIDHash.np)
	{
		fprintf(logfile,"call fill_PIDHash() before fresh_ID2Index()!\n");
		exit(1);
	}
	if(FRSH_GRPCAT==src_len)//cat
	{
		#ifndef GRPINPUT_INDEX //fofcat of JING's data are originally PIndex rather than PID
		Cat=(CATALOGUE *)src;
		for(i=0;i<Cat->Nids;i++)
			Cat->PIDorIndex[i]=lookup_ID2Ind(Cat->PIDorIndex[i]);//Now PIDorIndex has been refilled with particle Index in Pdat
		#endif	
	}
	else if(FRSH_SUBCAT==src_len)//subcat
	{
		SubCat=(SUBCATALOGUE *)src;
		/*====refresh SubCat with particle Index===*/
		for(i=0;i<SubCat->Nsubs;i++)
		{
			for(j=0;j<SubCat->SubLen[i];j++)
			{
				pid=SubCat->PSubArr[i][j];
				SubCat->PSubArr[i][j]=lookup_ID2Ind(pid);
			}
		}
	}
	else if(FRSH_SRCCAT==src_len)//srccat
	{
		SrcCat=(SRCCATALOGUE *)src;
		/*====refresh SubCat with particle Index===*/
		for(i=0;i<SrcCat->Nsubs;i++)
		{
			for(j=0;j<SrcCat->SubLen[i];j++)
			{
				pid=SrcCat->PSubArr[i][j];
				SrcCat->PSubArr[i][j]=lookup_ID2Ind(pid);
			}
			for(j=0;j<SrcCat->SubLen2[i];j++)
			{
				pid=SrcCat->PSubArr2[i][j];
				SrcCat->PSubArr2[i][j]=lookup_ID2Ind(pid);
			}
		}
	}
	else if(src_len>=0)
	{
		Arr=(HBTInt *)src;
		for(i=0;i<src_len;i++)
			Arr[i]=lookup_ID2Ind(Arr[i]);
	}
	else
	{
		 fprintf(logfile,"Error: wrong srclen when using fresh_ID2Index\nsrclen must be an integer no smaller than -3\n");fflush(logfile);
		exit(1);
	}
}

void create_sub_cat(SUBCATALOGUE *SubCat)//patched version
{
	HBTInt i,Nsubs;
	//SubCat->Ngroups=Ngroups;
	//SubCat->GrpOffset_Sub=mymalloc(sizeof(HBTInt)*Ngroups);
	//SubCat->GrpLen_Sub=mymalloc(sizeof(HBTInt)*Ngroups);//this block is added manually
	Nsubs=SubCat->Nsubs;	
	SubCat->SubLen=mymalloc(sizeof(HBTInt)*Nsubs);
	SubCat->SubOffset=mymalloc(sizeof(HBTInt)*Nsubs);
	SubCat->SubRank=mymalloc(sizeof(HBTInt)*Nsubs);
	SubCat->HaloChains=mymalloc(sizeof(struct Chain_data)*Nsubs);
	SubCat->sub_hierarchy=mymalloc(sizeof(struct Hierarchy)*Nsubs);
	SubCat->Property=mymalloc(sizeof(struct SubProperty)*Nsubs);
	SubCat->PSubArr=mymalloc(sizeof(HBTInt *)*Nsubs);
	for(i=0;i<Nsubs;i++)
		SubCat->PSubArr[i]=NULL;
}

void create_src_cat(SRCCATALOGUE *SrcCat)
{
	HBTInt i,Nsubs;
	Nsubs=SrcCat->Nsubs;//this must be assigned before using create()
	SrcCat->SubLen=mymalloc(sizeof(HBTInt)*Nsubs);
	SrcCat->SubLen2=mymalloc(sizeof(HBTInt)*Nsubs);
	SrcCat->SubOffset=mymalloc(sizeof(HBTInt)*Nsubs);
	SrcCat->CoreFrac=mymalloc(sizeof(HBTReal)*Nsubs);
	SrcCat->PSubArr=mymalloc(sizeof(HBTInt *)*Nsubs);
	SrcCat->PSubArr2=mymalloc(sizeof(HBTInt *)*Nsubs);
	for(i=0;i<Nsubs;i++)
	{
		SrcCat->PSubArr[i]=NULL;
		SrcCat->PSubArr2[i]=NULL;
	}
	SrcCat->HaloChains=NULL;
}

void save_sub_catalogue(HBTInt Nsnap, SUBCATALOGUE *Cat,char *SubCatPath)
{
FILE *fd;
char buf[1024];
HBTInt i,j;

  sprintf(buf, "%s/subcat_%03d", SubCatPath, (int)Nsnap);
  if(!(fd = fopen(buf, "w")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fwrite(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fwrite(Cat->GrpLen_Sub, sizeof(HBTInt), Cat->Ngroups, fd);
  fwrite(Cat->GrpOffset_Sub,sizeof(HBTInt), Cat->Ngroups, fd);
  fwrite(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fwrite(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
  fwrite(Cat->SubRank,sizeof(HBTInt), Cat->Nsubs, fd);
  fwrite(Cat->HaloChains,sizeof(struct Chain_data), Cat->Nsubs, fd);
  fwrite(Cat->sub_hierarchy,sizeof(struct Hierarchy),Cat->Nsubs,fd);
  fwrite(Cat->Property,sizeof(struct SubProperty), Cat->Nsubs, fd);

for(i=0;i<Cat->Nsubs;i++)
{
  #ifdef SHIFT_PID_0_1  //need to restore PID to starting from 0
  for(j=0;j<Cat->SubLen[i];j++)
		Cat->PSubArr[i][j]--;
  #endif
  fwrite(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
}
fwrite(&Cat->Nbirth,sizeof(HBTInt),1,fd);
fwrite(&Cat->NQuasi,sizeof(HBTInt),1,fd);
fwrite(&Cat->Ndeath,sizeof(HBTInt),1,fd);
fwrite(&Cat->Nsplitter,sizeof(HBTInt),1,fd);

  fclose(fd);
}

void load_sub_catalogue(HBTInt Nsnap, SUBCATALOGUE *Cat, char *SubCatPath)
{
FILE *fd;
char buf[1024];
HBTInt i,j;

  sprintf(buf, "%s/subcat_%03d", SubCatPath, (int)Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
    Cat->GrpOffset_Sub=mymalloc(sizeof(HBTInt)*Cat->Ngroups);   
	Cat->GrpLen_Sub=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  fread(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
	Cat->SubLen=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->SubOffset=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->SubRank=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->HaloChains=mymalloc(sizeof(struct Chain_data)*Cat->Nsubs);
	Cat->sub_hierarchy=mymalloc(sizeof(struct Hierarchy)*Cat->Nsubs);
	Cat->Property=mymalloc(sizeof(struct SubProperty)*Cat->Nsubs);
	Cat->PSubArr=mymalloc(sizeof(HBTInt *)*Cat->Nsubs);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fread(Cat->GrpLen_Sub, sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->GrpOffset_Sub,sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubRank,sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->HaloChains,sizeof(struct Chain_data), Cat->Nsubs, fd);
  fread(Cat->sub_hierarchy,sizeof(struct Hierarchy),Cat->Nsubs,fd);
  fread(Cat->Property,sizeof(struct SubProperty), Cat->Nsubs, fd);

for(i=0;i<Cat->Nsubs;i++)
{
	Cat->PSubArr[i]=mymalloc(sizeof(HBTInt)*Cat->SubLen[i]);
	fread(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
  #ifdef SHIFT_PID_0_1  //now shift PID from 0 to 1
  for(j=0;j<Cat->SubLen[i];j++)
		Cat->PSubArr[i][j]++;
  #endif
}
fread(&Cat->Nbirth,sizeof(HBTInt),1,fd);
fread(&Cat->NQuasi,sizeof(HBTInt),1,fd);
fread(&Cat->Ndeath,sizeof(HBTInt),1,fd);
fread(&Cat->Nsplitter,sizeof(HBTInt),1,fd);

  fclose(fd);
}
void load_sub_table(HBTInt Nsnap, SUBCATALOGUE *Cat, char *SubCatPath)
{//only load global properties, without filling particles
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/subcat_%03d", SubCatPath, (int)Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
    Cat->GrpOffset_Sub=mymalloc(sizeof(HBTInt)*Cat->Ngroups);   
	Cat->GrpLen_Sub=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  fread(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
	Cat->SubLen=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->SubOffset=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->SubRank=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->HaloChains=mymalloc(sizeof(struct Chain_data)*Cat->Nsubs);
	Cat->sub_hierarchy=mymalloc(sizeof(struct Hierarchy)*Cat->Nsubs);
	Cat->Property=mymalloc(sizeof(struct SubProperty)*Cat->Nsubs);
	Cat->PSubArr=NULL;
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fread(Cat->GrpLen_Sub, sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->GrpOffset_Sub,sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubRank,sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->HaloChains,sizeof(struct Chain_data), Cat->Nsubs, fd);
  fread(Cat->sub_hierarchy,sizeof(struct Hierarchy),Cat->Nsubs,fd);
  fread(Cat->Property,sizeof(struct SubProperty), Cat->Nsubs, fd);

	fseek(fd,sizeof(HBTInt)*Cat->Nids,SEEK_CUR);
//~ for(i=0;i<Cat->Nsubs;i++)
//~ {
	//~ Cat->PSubArr[i]=mymalloc(sizeof(HBTInt)*Cat->SubLen[i]);
	//~ fread(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
//~ }
fread(&Cat->Nbirth,sizeof(HBTInt),1,fd);
fread(&Cat->NQuasi,sizeof(HBTInt),1,fd);
fread(&Cat->Ndeath,sizeof(HBTInt),1,fd);
fread(&Cat->Nsplitter,sizeof(HBTInt),1,fd);

  fclose(fd);
}
void free_sub_table(SUBCATALOGUE *SubCat)//corresponding to load_sub_table
{
	myfree(SubCat->GrpOffset_Sub);
	myfree(SubCat->GrpLen_Sub);
	myfree(SubCat->SubLen);
	myfree(SubCat->SubOffset);
	myfree(SubCat->SubRank);
	myfree(SubCat->HaloChains);
	myfree(SubCat->sub_hierarchy);
	myfree(SubCat->Property);
}
void save_src_catalogue(HBTInt Nsnap, SRCCATALOGUE *Cat,char *SrcCatPath)
{
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/srccat_%03d", SrcCatPath, (int)Nsnap);
  if(!(fd = fopen(buf, "w")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fwrite(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fwrite(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fwrite(Cat->SubLen2, sizeof(HBTInt), Cat->Nsubs, fd);
  fwrite(Cat->CoreFrac, sizeof(HBTReal), Cat->Nsubs, fd);
  fwrite(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
for(i=0;i<Cat->Nsubs;i++)
  fwrite(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
for(i=0;i<Cat->Nsubs;i++)
{
	if(Cat->SubLen2[i]>0)
 	 fwrite(Cat->PSubArr2[i], sizeof(HBTInt), Cat->SubLen2[i], fd);
}
fwrite(&Cat->NDeathSp,sizeof(HBTInt),1,fd);
  fclose(fd);
}
void load_src_catalogue(HBTInt Nsnap, SRCCATALOGUE *Cat,char *SrcCatPath)
{
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/srccat_%03d", SrcCatPath, (int)Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  Cat->SubLen=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
  Cat->SubLen2=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
  Cat->CoreFrac=mymalloc(sizeof(HBTReal)*Cat->Nsubs);
  Cat->SubOffset=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
  fread(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubLen2, sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->CoreFrac, sizeof(HBTReal), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
  Cat->PSubArr=mymalloc(sizeof(HBTInt*)*Cat->Nsubs);
for(i=0;i<Cat->Nsubs;i++)
{
	Cat->PSubArr[i]=mymalloc(sizeof(HBTInt)*Cat->SubLen[i]);
  	fread(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
}
  Cat->PSubArr2=mymalloc(sizeof(HBTInt*)*Cat->Nsubs);
for(i=0;i<Cat->Nsubs;i++)
{
	if(Cat->SubLen2[i]>0)
	{
	Cat->PSubArr2[i]=mymalloc(sizeof(HBTInt)*Cat->SubLen2[i]);
  	fread(Cat->PSubArr2[i], sizeof(HBTInt), Cat->SubLen2[i], fd);
	}
	else
	Cat->PSubArr2[i]=NULL;
}
fread(&Cat->NDeathSp,sizeof(HBTInt),1,fd);
  fclose(fd);
  Cat->HaloChains=NULL;
}

void load_sp2pro(HBTInt Nsnap_dest,HBTInt *P2Npro,HBTInt *P2Nsplitter_dest, HBTInt **P2sp2pro,char *SubcatPath)
{/*	Npro=SubCatPro.Nsubs;	Nsplitter_dest=SubCatDest.Nsplitter
	* sp2pro(ProSubID) fixes ProSubID to the valid subid if ProSubID>=Npro;
	* returns NULL if Nsplitter==0
	* 
here you have two adjacent snapshots, say snap_pro and snap_dest.
there would be Npro subhalos in snap_pro, some of which may split and some not. 
* but the total number of splitted out pieces (not including the main part--let's name it main splinter--which is marked as the official descendant of this progenitor subhalo) are Nsplitters. 
* These Nsplitter pieces only get to be identified as subhalos in snap_dest, and they share physical progenitors with the main splinters. 
* To be able to tell that they are just splitted pieces from their progenitors, 
* we do not record their progenitor's physical subhaloID as their ProSubID,
*  but create some virtual subhaloIDs as their ProSubID.  
* These virtual ProSubIDs are in the range [Npro, Npro+Nsplitters-1] ,
*  so you can NOT directly find out an existing subhalo from snap_pro as these splinters' progenitor.

Then sp2pro comes to help you translate a virtual ProSubID to the physical SubID which creates the splitter. sp2pro just stores the physical progenitor SubID for each virtual ProSubID. You see there would be Nsplitters elements in this array. In fortran you may want to create this array with subscript ranging from Npro to Npro+Nsplitters-1, so that you can get real subids by calling
RealSubID=sp2pro[ProSubID]
In C, this is similarly achieved by moving the pointer as
sp2pro=sp2pro-Npro,
so that you can use RealSubID=sp2pro(ProSubID).

P2***  are pointers to variables. By using pointers, or address of the input variables, we can directly modify the input variables. So in fact these P2*** are all used for outputing these variables.
	* */
	
	char buf[1024];
	HBTInt *sp2pro,Npro,Nsplitter;
	FILE *fp;
	
	sprintf(buf, "%s/splitters/sp2pro_%03d",SubcatPath, (int)Nsnap_dest);
	if(!(fp = fopen(buf, "r")))
	{
		fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
		exit(1);
	}
	fread(&Nsplitter,sizeof(HBTInt),1,fp);
	fread(&Npro,sizeof(HBTInt),1,fp);
	*P2Npro=Npro;
	*P2Nsplitter_dest=Nsplitter;
	*P2sp2pro=NULL;
	if(Nsplitter)
	{
	sp2pro=mymalloc(sizeof(HBTInt)*Nsplitter);
	fread(sp2pro,sizeof(HBTInt),Nsplitter,fp);
	*P2sp2pro=sp2pro-Npro;
	}
	fclose(fp);
}
void free_sp2pro(HBTInt *sp2pro,HBTInt Npro,HBTInt Nsplitter_dest)
{
	if(Nsplitter_dest)
	free(sp2pro+Npro);
}
void free_src_catalogue(SRCCATALOGUE *SrcCat)
{
	myfree(SrcCat->SubLen);
	myfree(SrcCat->SubLen2);
	myfree(SrcCat->CoreFrac);
	myfree(SrcCat->SubOffset);
	myfree(SrcCat->PSubArr);
	myfree(SrcCat->PSubArr2);
	myfree(SrcCat->HaloChains);
}
void free_sub_catalogue(SUBCATALOGUE *SubCat)//patched version
{
	myfree(SubCat->GrpOffset_Sub);
	myfree(SubCat->GrpLen_Sub);
	myfree(SubCat->SubLen);
	myfree(SubCat->SubOffset);
	myfree(SubCat->SubRank);
	myfree(SubCat->HaloChains);
	myfree(SubCat->sub_hierarchy);
	myfree(SubCat->Property);
	myfree(SubCat->PSubArr);//this only frees the pointer array, remember to reuse or free the actual pid array.
}
void erase_sub_catalogue(SUBCATALOGUE *SubCat)
{
	HBTInt i;
	for(i=0;i<SubCat->Nsubs;i++)
	{
		myfree(SubCat->PSubArr[i]);
	}
	free_sub_catalogue(SubCat);
}
void erase_src_catalogue(SRCCATALOGUE *SrcCat)
{
	HBTInt i;
	for(i=0;i<SrcCat->Nsubs;i++)
	{
		myfree(SrcCat->PSubArr[i]);
		myfree(SrcCat->PSubArr2[i]);
	}
	free_src_catalogue(SrcCat);	
}
void fix_spfile(FILE **p2fp,char *buf,HBTInt Snap)
{
	FILE *fp1,*fp2;
	char buf2[1024],cmd[1024];
	int tmp[4],nread,line;
	fp1=*p2fp;
	sprintf(buf2,"%s._tmp",buf);
	myfopen(fp2,buf2,"w");
	fprintf(fp1,"%d,"HBTIFMT",%d,%d,-1\n",-1,Snap,-1,-1);//add safe end
	fclose(fp1);
	myfopen(fp1,buf,"r");//since fpsp was opened with "a" mode which can not be read, reopen it with "r" here.
	fscanf(fp1,"%d,%d,%d,%d\n",tmp,tmp+1,tmp+2,tmp+3);
	//~ printf("%d,%d,%d,%d\n",tmp[0],tmp[1],tmp[2],tmp[3]);
	line=1;
	while(tmp[1]<Snap-1)//[snap-1,snap] interval  will be generated for current loop,so eliminate any previous remnent from snap-1
	{
	do{
			fprintf(fp2,"%d,%d,%d,%d\n",tmp[0],tmp[1],tmp[2],tmp[3]);	
			nread=fscanf(fp1,"%d,%d,%d,%d\n",tmp,tmp+1,tmp+2,tmp+3);
				if(nread!=4||ferror(fp1))
				{
					fprintf(logfile,"error reading splitter file %s,at line %d,%d\n",buf,line+1,nread);fflush(logfile);
					exit(1);
				}
			line++;
		}
	while(tmp[0]>=0);
	}
	fclose(fp1);
	fclose(fp2);
	sprintf(cmd,"rm %s; mv %s %s\n",buf,buf2,buf);
	system(cmd);
	myfopen(*p2fp,buf,"a");
}

void load_pro2dest(HBTInt Nsnap_pro,HBTInt **P2pro2dest,HBTInt *P2Nsubs,char *destdir)
{
FILE *fp;
char buf[1024];
HBTInt Nsubs,*pro2dest;

  sprintf(buf, "%s/pro2dest/pro2dest_%03d", destdir, (int)Nsnap_pro);
  myfopen(fp,buf,"r");
  fread(&Nsubs,sizeof(HBTInt),1,fp);//Nsubs=Npro+Nspl
  pro2dest=mymalloc(sizeof(HBTInt)*(Nsubs+1));
  pro2dest++;
  fread(pro2dest,sizeof(HBTInt),Nsubs,fp);
  fclose(fp);
  pro2dest[-1]=-1;//no pro, no dest
  *P2Nsubs=Nsubs;
  *P2pro2dest=pro2dest;
}
void free_pro2dest(HBTInt *pro2dest)
{
	free(pro2dest-1);
}
void save_pro2dest(HBTInt Nsnap_pro,HBTInt *pro2dest,HBTInt Nsubs,char *destdir)
{
FILE *fp;
char buf[1024];

  sprintf(buf, "%s/pro2dest/pro2dest_%03d", destdir, (int)Nsnap_pro);
  myfopen(fp,buf,"w");
  fwrite(&Nsubs,sizeof(HBTInt),1,fp);
  fwrite(pro2dest,sizeof(HBTInt),Nsubs,fp);
  fclose(fp);
}	

void save_group_catalogue_HBT(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{//HBT groupcat format, same as old gadget groups, except using HBTInt datatype 
  FILE *fd;
  HBTInt i,NFiles;
  char buf[1024];
  
#ifdef SNAPLIST
Nsnap=snaplist[Nsnap];
#endif
  sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
  myfopen(fd,buf,"w");

  fwrite(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  NFiles=1;
  fwrite(&NFiles, sizeof(HBTInt), 1, fd);
  
  fwrite(Cat->Len, sizeof(HBTInt), Cat->Ngroups, fd);
  fwrite(Cat->Offset,sizeof(HBTInt), Cat->Ngroups, fd);
  fclose(fd);

  sprintf(buf, "%s/group_ids_%03d", GrpPath, (int)Nsnap);
  myfopen(fd,buf,"w");

  fwrite(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fwrite(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  NFiles=1;
  fwrite(&NFiles, sizeof(HBTInt), 1, fd);

  HBTInt *PIDs;
  PIDs=mymalloc(sizeof(HBTInt)*Cat->Nids);
  for(i=0;i<Cat->Nids;i++) 
	PIDs[i]=lookup_Ind2ID(Cat->PIDorIndex[i]);
  fwrite(PIDs,sizeof(HBTInt),Cat->Nids,fd);
  myfree(PIDs);
  
  fclose(fd);
}

void load_group_catalogue_HBT(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{//old groupcat format
  FILE *fd;
  HBTInt i,NFiles;
  char buf[1024];
  
#ifdef SNAPLIST
Nsnap=snaplist[Nsnap];
#endif
  sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&NFiles, sizeof(HBTInt), 1, fd);
  
  Cat->Len=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  fread(Cat->Len, sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->Offset,sizeof(HBTInt), Cat->Ngroups, fd);
  fclose(fd);

fprintf(logfile,"Snap="HBTIFMT" Ngroups="HBTIFMT" Nids="HBTIFMT"\n", Nsnap, Cat->Ngroups, Cat->Nids);
  sprintf(buf, "%s/group_ids_%03d", GrpPath, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&NFiles, sizeof(HBTInt), 1, fd);

  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
  fread(Cat->PIDorIndex,sizeof(HBTInt),Cat->Nids,fd);
    
  fclose(fd);
  
  Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);//consider move this out.............
}

/*==Bug fixes/Changes==*
 * bug: fix_spfile() Snap duplication 25.03.2009*
 * bug: fresh_ID2Index missed to fresh SrcCat.PSubArr2[][] 15.05.2009
 * ======*/
	
