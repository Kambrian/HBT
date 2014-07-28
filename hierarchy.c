#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

void transfer_subcat(SUBCATALOGUE *SubCatTo,SUBCATALOGUE *SubCatFrom)
{
	free_sub_catalogue(SubCatTo);
	memcpy(SubCatTo,SubCatFrom,sizeof(SUBCATALOGUE));
	//Nullify SubCatFrom (buf not freeing the storage since they have been taken over by SubTo)
	SubCatFrom->Ngroups=0;
	SubCatFrom->Nsubs=0;
	SubCatFrom->Nids=0;
	SubCatFrom->GrpLen_Sub=NULL;
	SubCatFrom->GrpOffset_Sub=NULL;
	SubCatFrom->SubLen=NULL;
	SubCatFrom->SubOffset=NULL;
	SubCatFrom->SubRank=NULL;
	SubCatFrom->HaloChains=NULL;
	SubCatFrom->Property=NULL;
	SubCatFrom->sub_hierarchy=NULL;
	SubCatFrom->PSubArr=NULL;
	SubCatFrom->Nbirth=0;
	SubCatFrom->NQuasi=0;
	SubCatFrom->Ndeath=0;
	SubCatFrom->Nsplitter=0;
}
void transfer_srccat(SRCCATALOGUE *SrcCatTo,SRCCATALOGUE *SrcCatFrom)
{
	free_src_catalogue(SrcCatTo);
	memcpy(SrcCatTo,SrcCatFrom,sizeof(SRCCATALOGUE));
	//Nullify SrcCatFrom
	SrcCatFrom->Nsubs=0;
	SrcCatFrom->Nids=0;
	SrcCatFrom->SubLen=NULL;
	SrcCatFrom->SubLen2=NULL;
	SrcCatFrom->CoreFrac=NULL;
	SrcCatFrom->SubOffset=NULL;
	SrcCatFrom->HaloChains=NULL;
	SrcCatFrom->PSubArr=NULL;
	SrcCatFrom->PSubArr2=NULL;
	SrcCatFrom->NDeathSp=0;
}

void break_out_sub_family(HBTInt mainsubID,SUBCATALOGUE *SubCat)
{
	HBTInt sib,nibs,pre;
	if((nibs=SubCat->sub_hierarchy[mainsubID].nibs)>=0)//only mainsubs have nibs<0; if a main, no need to breakout
	{
		pre=SubCat->sub_hierarchy[mainsubID].pre;
		sib=SubCat->sub_hierarchy[mainsubID].next;
		if(pre>=0)
		{
			SubCat->sub_hierarchy[pre].next=sib;
			SubCat->sub_hierarchy[mainsubID].pre=-1;
		}
		else 
			SubCat->sub_hierarchy[nibs].sub=sib;
		if(sib>=0)
		{
			SubCat->sub_hierarchy[sib].pre=pre;
			SubCat->sub_hierarchy[mainsubID].next=-1;
		}
		SubCat->sub_hierarchy[mainsubID].nibs=-1;
	}
}

void markcross_N_kickdeath_recursive(HBTInt mainsubID,SUBCATALOGUE *SubCat)
{
	/*SubCat->HaloChains has to be raw,that is, sorted so that its subscript=proID*/
	HBTInt son,sib,nibs,pre;	

	if(SubCat->SubLen[mainsubID]==0)//this sub dies
	{
		if((son=SubCat->sub_hierarchy[mainsubID].sub)>=0)//level-up its sons and delete itself
		{	
			SubCat->sub_hierarchy[mainsubID].sub=-1;
			if((nibs=SubCat->sub_hierarchy[mainsubID].nibs)<0)//a main,break-up all its sons 
			{
			do{
					sib=son;
					SubCat->sub_hierarchy[sib].nibs=-1;
					SubCat->sub_hierarchy[sib].pre=-1;
					son=SubCat->sub_hierarchy[sib].next;
					SubCat->sub_hierarchy[sib].next=-1;
					markcross_N_kickdeath_recursive(sib,SubCat);//and mark_subcross
				}while(son>=0);
			}
			else//need left or upper and perhaps right link
			{	
				SubCat->sub_hierarchy[mainsubID].nibs=-1;	
				SubCat->sub_hierarchy[son].nibs=nibs;
				if((pre=SubCat->sub_hierarchy[mainsubID].pre)>0)//link left
				{
				SubCat->sub_hierarchy[mainsubID].pre=-1;
				SubCat->sub_hierarchy[son].pre=pre;
				SubCat->sub_hierarchy[pre].next=son;
				}
				else// link up
				{
					SubCat->sub_hierarchy[nibs].sub=son;
				}
				sib=SubCat->sub_hierarchy[son].next;
				while(sib>=0)// link others up
				{
					markcross_N_kickdeath_recursive(son,SubCat);//update last sub
					son=sib;
					SubCat->sub_hierarchy[son].nibs=nibs;
					sib=SubCat->sub_hierarchy[son].next;
				}
				if((sib=SubCat->sub_hierarchy[mainsubID].next)>=0)//link right,welcome new sib
				{
					SubCat->sub_hierarchy[mainsubID].next=-1;
					SubCat->sub_hierarchy[son].next=sib;
					SubCat->sub_hierarchy[sib].pre=son;
				}
				markcross_N_kickdeath_recursive(son,SubCat);//update last sub
			}
		}
		else break_out_sub_family(mainsubID,SubCat);//no son,but have siblings to link
	}
	else 	if((son=SubCat->sub_hierarchy[mainsubID].sub)>=0)//alive, so update its sons
	{	
		sib=SubCat->sub_hierarchy[son].next;
		if(SubCat->HaloChains[mainsubID].HostID!=SubCat->HaloChains[son].HostID) //level-up this sub if it crossed-out
		{
			SubCat->sub_hierarchy[mainsubID].sub=sib; // no pre but nibs to link
			if(sib>=0)
			{
				SubCat->sub_hierarchy[sib].pre=-1;	
				markcross_N_kickdeath_recursive(mainsubID,SubCat);
			}
			SubCat->sub_hierarchy[son].nibs=-1;
			SubCat->sub_hierarchy[son].next=-1;
			
			markcross_N_kickdeath_recursive(son,SubCat);		
		}
		else //the big son stayed, so its sibs got "pre" to link
		{
			markcross_N_kickdeath_recursive(son,SubCat);//update big son
			while(sib>=0)
			{
				son=sib;
				sib=SubCat->sub_hierarchy[son].next;
				if(SubCat->HaloChains[mainsubID].HostID!=SubCat->HaloChains[son].HostID)	//level-up this sub if it crossed-out
				{
					SubCat->sub_hierarchy[SubCat->sub_hierarchy[son].pre].next=sib;
					if(sib>=0)
					{
						SubCat->sub_hierarchy[sib].pre=SubCat->sub_hierarchy[son].pre;
					}
					SubCat->sub_hierarchy[son].nibs=-1;
					SubCat->sub_hierarchy[son].next=-1;
					SubCat->sub_hierarchy[son].pre=-1;
				}
				markcross_N_kickdeath_recursive(son,SubCat);
			}
		}
	}
}
void kickdeath_recursive(HBTInt mainsubID,SUBCATALOGUE *SubCat)
{	/*simplified version of markdeath_N_kickdeath_recursive(), to only kick death*/
	/*SubCat->HaloChains has to be raw,that is, sorted so that its subscript=proID*/
	HBTInt son,sib,nibs,pre;	

	if(SubCat->SubLen[mainsubID]==0)//this sub dies
	{
		if((son=SubCat->sub_hierarchy[mainsubID].sub)>=0)//level-up its sons and delete itself
		{	
			SubCat->sub_hierarchy[mainsubID].sub=-1;
			if((nibs=SubCat->sub_hierarchy[mainsubID].nibs)<0)//a main,break-up all its sons 
			{
			do{
					sib=son;
					SubCat->sub_hierarchy[sib].nibs=-1;
					SubCat->sub_hierarchy[sib].pre=-1;
					son=SubCat->sub_hierarchy[sib].next;
					SubCat->sub_hierarchy[sib].next=-1;
					kickdeath_recursive(sib,SubCat);//and mark_subcross
				}while(son>=0);
			}
			else//need left or upper and perhaps right link
			{	
				SubCat->sub_hierarchy[mainsubID].nibs=-1;	
				SubCat->sub_hierarchy[son].nibs=nibs;
				if((pre=SubCat->sub_hierarchy[mainsubID].pre)>0)//link left
				{
				SubCat->sub_hierarchy[mainsubID].pre=-1;
				SubCat->sub_hierarchy[son].pre=pre;
				SubCat->sub_hierarchy[pre].next=son;
				}
				else// link up
				{
					SubCat->sub_hierarchy[nibs].sub=son;
				}
				sib=SubCat->sub_hierarchy[son].next;
				while(sib>=0)// link others up
				{
					kickdeath_recursive(son,SubCat);//update last sub
					son=sib;
					SubCat->sub_hierarchy[son].nibs=nibs;
					sib=SubCat->sub_hierarchy[son].next;
				}
				if((sib=SubCat->sub_hierarchy[mainsubID].next)>=0)//link right,welcome new sib
				{
					SubCat->sub_hierarchy[mainsubID].next=-1;
					SubCat->sub_hierarchy[son].next=sib;
					SubCat->sub_hierarchy[sib].pre=son;
				}
				kickdeath_recursive(son,SubCat);//update last sub
			}
		}
		else break_out_sub_family(mainsubID,SubCat);//no son,but have siblings to link
	}
	else//alive,still need to check its sons
	{
		if((son=SubCat->sub_hierarchy[mainsubID].sub)>=0)//kick sons if any
		{	
			do{
					sib=son;
					son=SubCat->sub_hierarchy[sib].next;
					kickdeath_recursive(sib,SubCat);//and mark_subcross
				}while(son>=0);
		}	
	}
}

static int compare_splitterlen(const void *a, const void *b)//used to sort desendent id in descending order
{
  if(((struct cand_data *) a)->SubLen > ((struct cand_data *) b)->SubLen)
    return -1;

  if(((struct cand_data *) a)->SubLen < ((struct cand_data *) b)->SubLen)
    return +1;

  return 0;
}
static HBTInt *ID2HaloCMP;
static int compare_id2halo(const void *a, const void *b)//used to sort particle id in ascending order of it's host
{
  if(ID2HaloCMP[*(HBTInt *)a] > ID2HaloCMP[*(HBTInt *)b])
    return +1;

  if(ID2HaloCMP[*(HBTInt *)a] < ID2HaloCMP[*(HBTInt *)b])
    return -1;

  return 0;
}
void PARAsplit_srccat(CATALOGUE* Cat,SRCCATALOGUE *SrcCat, HBTInt SnapshotNum)
{	//parallel version of split_srccat(),using qsort;
	/*split SrcCat, and make HaloChains
	 *        For MAJOR split (main splitter fraction is smaller than 0.7),
	 *            those splitted subs will be reset of their SubLen2=-1 and CoreFrac=CoreFrac0 and PSubArr2=NULL
	 *             this is also the case when sublen2>sublen_new/MassRelax,as if sublen2 haven't been picked up
	 * 		otherwise sublen2 and PSubArr2 and CoreFrac keep intact
	 * Cat->ID2halo will be filled 
	 * */
	//these automatic vars are all private implicitly!!!
	FILE *fpsp;
	char buf[1024];
	HBTInt i,splid,pid,proID,desID,Ncands,candNids,candOffset,maxcands,Nsubs_old;
	struct cand_data *SubSplArr;
	//declared as static here to make them shared vars during omp parallelization
	static HBTInt Nsplitter,Ndeathsp;
	static struct SubSpl
	{
		struct cand_data *SubSplArr;//SubSplArr[Ncands]
		HBTInt Ncands;
	}*CatSplArr;

	#ifdef HALO_PARA
	#pragma omp single
	#endif
	{
	ID2HaloCMP=Cat->ID2Halo;
	Nsplitter=0;Ndeathsp=0;
	CatSplArr=mymalloc(sizeof(struct SubSpl)*SrcCat->Nsubs);
	SrcCat->HaloChains=mymalloc(sizeof(struct Chain_data)*SrcCat->Nsubs);//HaloChains should be freed or passed away before allocation
	}
	#ifdef HALO_PARA
	#pragma omp for reduction(+: Ndeathsp,Nsplitter) schedule(dynamic,1)
	#endif
	for(proID=0;proID<SrcCat->Nsubs;proID++)//find hosts for the bound structure
	{	
		//~ HaloLinkCount=mymalloc(sizeof(HBTInt)*(Cat->Ngroups+1));//will contain the backgroud as an extra halo
		//~ HaloLinkCount++;
		SrcCat->HaloChains[proID].ProSubID=proID;
		qsort(SrcCat->PSubArr[proID],SrcCat->SubLen[proID],sizeof(HBTInt),compare_id2halo);
		maxcands=20;
		CatSplArr[proID].SubSplArr=mymalloc(sizeof(struct cand_data)*maxcands);
		desID=-2;candNids=0;candOffset=0;Ncands=0;
		for(i=0;i<SrcCat->SubLen[proID];i++)		//find desendents:
		{
			pid=SrcCat->PSubArr[proID][i];
			if(desID!=ID2HaloCMP[pid])//found new host
			{
				if(candNids>=NSRCMIN)//previous cand is a splitter
				{
					CatSplArr[proID].SubSplArr[Ncands].desID=desID;
					CatSplArr[proID].SubSplArr[Ncands].SubLen=candNids;
					CatSplArr[proID].SubSplArr[Ncands].SubArr=mymalloc(sizeof(HBTInt)*candNids);
					memcpy(CatSplArr[proID].SubSplArr[Ncands].SubArr,SrcCat->PSubArr[proID]+candOffset,sizeof(HBTInt)*candNids);
					Ncands++;
					if(Ncands==maxcands)
					{
						fprintf(logfile,"maxcands>"HBTIFMT",proID="HBTIFMT",proLen="HBTIFMT"\n",maxcands,proID,SrcCat->SubLen[proID]);fflush(logfile);
						maxcands*=2;
						CatSplArr[proID].SubSplArr=realloc(CatSplArr[proID].SubSplArr,sizeof(struct cand_data)*maxcands);
					}
				}
			candOffset+=candNids;	
			desID=ID2HaloCMP[pid];
			candNids=1;
			}
			else
			candNids++;
		}
		if(candNids>=NSRCMIN)//previous cand is a splitter
		{
			CatSplArr[proID].SubSplArr[Ncands].desID=desID;
			CatSplArr[proID].SubSplArr[Ncands].SubLen=candNids;
			CatSplArr[proID].SubSplArr[Ncands].SubArr=mymalloc(sizeof(HBTInt)*candNids);
			memcpy(CatSplArr[proID].SubSplArr[Ncands].SubArr,SrcCat->PSubArr[proID]+candOffset,sizeof(HBTInt)*candNids);
			Ncands++;
		}
		CatSplArr[proID].Ncands=Ncands;
		
		if(Ncands>0)//split this sub
		{
		SubSplArr=CatSplArr[proID].SubSplArr;
		if(Ncands>1) qsort(SubSplArr,Ncands,sizeof(struct cand_data),compare_splitterlen);
		//move main splitters back to srccat
		SrcCat->HaloChains[proID].HostID=SubSplArr[0].desID;
		if((SubSplArr[0].SubLen<0.7*SrcCat->SubLen[proID])||SrcCat->SubLen2[proID]*MassRelax_Input>SubSplArr[0].SubLen)//splitted part is larger than 0.3, defined as a major split; or when Len2 is too big
		{//discard sublen2 and reset CoreFrac
			SrcCat->SubLen2[proID]=-1;
			myfree(SrcCat->PSubArr2[proID]);SrcCat->PSubArr2[proID]=NULL;
			SrcCat->CoreFrac[proID]=CoreFrac0;
		}			
		SrcCat->SubLen[proID]=SubSplArr[0].SubLen;
		free(SrcCat->PSubArr[proID]);
		SrcCat->PSubArr[proID]=SubSplArr[0].SubArr; 
		Nsplitter+=(Ncands-1);	
		}
		else
		{
			SrcCat->HaloChains[proID].HostID=-2;//splitted to death
			SrcCat->SubLen[proID]=0;//note: this would make problem when making trees; OK, solved.
			SrcCat->SubLen2[proID]=0;
			SrcCat->CoreFrac[proID]=CoreFrac0;
			free(SrcCat->PSubArr[proID]);
			SrcCat->PSubArr[proID]=NULL;
			if(SrcCat->PSubArr2[proID]!=NULL)	
			{
				free(SrcCat->PSubArr2[proID]);
				SrcCat->PSubArr2[proID]=NULL;
			}
			Ndeathsp++;
		}
	}
		#ifdef HALO_PARA
		#pragma omp single
		#endif
		{
		//~ printf("splitters: "HBTIFMT"\n",Nsplitter);
		//fill up splitters into subcattmp and halochains
		if(0==SrcCat->Nsubs)
		{
			Ndeathsp=0;
			Nsplitter=0;
		}//in case the un-executed reduction clause make Ndeathsp and Nsplitter uninitialized 
		SrcCat->NDeathSp=Ndeathsp;
		sprintf(buf, "%s/splitters/sp2pro_%03d",SUBCAT_DIR,(int)SnapshotNum);
		myfopen(fpsp,buf,"w");
		fwrite(&Nsplitter,sizeof(HBTInt),1,fpsp);
		fwrite(&SrcCat->Nsubs,sizeof(HBTInt),1,fpsp);
		Nsubs_old=SrcCat->Nsubs; 
		if(Nsplitter)
		{
		splid=SrcCat->Nsubs;
		SrcCat->Nsubs+=Nsplitter;
		SrcCat->PSubArr=realloc(SrcCat->PSubArr,sizeof(HBTInt *)*SrcCat->Nsubs);
		SrcCat->PSubArr2=realloc(SrcCat->PSubArr2,sizeof(HBTInt *)*SrcCat->Nsubs);
		SrcCat->SubLen=realloc(SrcCat->SubLen,sizeof(HBTInt)*SrcCat->Nsubs);
		SrcCat->SubLen2=realloc(SrcCat->SubLen2,sizeof(HBTInt)*SrcCat->Nsubs);
		SrcCat->CoreFrac=realloc(SrcCat->CoreFrac,sizeof(HBTReal)*SrcCat->Nsubs);
		SrcCat->HaloChains=realloc(SrcCat->HaloChains,sizeof(struct Chain_data)*SrcCat->Nsubs);
		for(proID=0;proID<Nsubs_old;proID++)
		{
			if(CatSplArr[proID].Ncands>1)
			{
				for(i=1;i<CatSplArr[proID].Ncands;i++)
				{
					SrcCat->PSubArr[splid]=CatSplArr[proID].SubSplArr[i].SubArr;
					SrcCat->SubLen[splid]=CatSplArr[proID].SubSplArr[i].SubLen;
					SrcCat->HaloChains[splid].HostID=CatSplArr[proID].SubSplArr[i].desID;
					SrcCat->HaloChains[splid].ProSubID=splid;
					SrcCat->PSubArr2[splid]=NULL;
					SrcCat->SubLen2[splid]=-1;
					SrcCat->CoreFrac[splid]=CoreFrac0;
					fwrite(&proID,sizeof(HBTInt),1,fpsp);
					splid++;
				}
			}
		}
		if(splid!=SrcCat->Nsubs){ fprintf(logfile,"error:splid and Nsubs mismatch:"HBTIFMT","HBTIFMT"\n",splid,SrcCat->Nsubs);exit(1);}
		}
		fclose(fpsp);
		for(proID=0;proID<Nsubs_old;proID++)
		free(CatSplArr[proID].SubSplArr);
		free(CatSplArr);
		}
}

void PARAmake_srcsub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB,SRCCATALOGUE *SrcCatA,SRCCATALOGUE *SrcCatB)
{/*  *output: SubCatA,SrcCatA
		*useful part: 
		* 						SubLen,SubArr, HaloChains from SrcCat
		* 						sub_hierarchy from SubCat and Nsplitter
		* 						expand SubCatA.CoM
		* 						Nsplitter from SrcCat-SubCat
		* 						Nsubs from SrcCat						                *
		* when returned, SubCatB and SrcCatB 's pointers needn't/mustn't be freed before allocation,their memory has been taken over by A*/
	HBTInt i;

#pragma omp single
{
	transfer_subcat(SubCatA,SubCatB);//free A, transfer B's data to A and nullify B
	SubCatA->HaloChains=SrcCatB->HaloChains; SrcCatB->HaloChains=NULL;//Take over the HaloChains here,and Nullify B to protect it,now only SubA has a chain
	SubCatA->Nsplitter=SrcCatB->Nsubs-SubCatA->Nsubs;
}
	if(SubCatA->Nsplitter)
	{	
	#pragma omp single	
		SubCatA->sub_hierarchy=realloc(SubCatA->sub_hierarchy,sizeof(struct Hierarchy)*SrcCatB->Nsubs);
	#pragma omp for 
		for(i=SubCatA->Nsubs;i<SrcCatB->Nsubs;i++)//fill splitter.sub_hierarchy
		{
			SubCatA->sub_hierarchy[i].nibs=-1;
			SubCatA->sub_hierarchy[i].pre=-1;
			SubCatA->sub_hierarchy[i].next=-1;
			SubCatA->sub_hierarchy[i].sub=-1;
		}
	#pragma omp single
	{
		free(SubCatA->SubLen);
		SubCatA->SubLen=mymalloc(sizeof(HBTInt)*SrcCatB->Nsubs);
		free(SubCatA->Property);
		SubCatA->Property=mymalloc(sizeof(struct SubProperty)*SrcCatB->Nsubs);
		SubCatA->Nsubs=SrcCatB->Nsubs;
	}
	}
	#pragma omp single
	memcpy(SubCatA->SubLen,SrcCatB->SubLen,sizeof(HBTInt)*SrcCatB->Nsubs);//copy the arrs
	#pragma omp for 
	for(i=0;i<SubCatA->Nsubs-SubCatA->Nsplitter;i++)
		free(SubCatA->PSubArr[i]);
	#pragma omp single
	{
	myfree(SubCatA->PSubArr);
	SubCatA->PSubArr=mymalloc(sizeof(HBTInt *)*SrcCatB->Nsubs);
	}
	#pragma omp for
	for(i=0;i<SrcCatB->Nsubs;i++)//copy the ids
	{
		SubCatA->PSubArr[i]=mymalloc(sizeof(HBTInt)*SrcCatB->SubLen[i]);
		memcpy(SubCatA->PSubArr[i],SrcCatB->PSubArr[i],sizeof(HBTInt)*SrcCatB->SubLen[i]);
	}
	#pragma omp single
	transfer_srccat(SrcCatA,SrcCatB);
}

void migrate_sub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt SubRank,HBTInt HostID,HBTInt *pro2dest)
{
	HBTInt desSubID;
	desSubID=pro2dest[proSubID];
	SubCatB->PSubArr[desSubID]=SubCatA->PSubArr[proSubID];SubCatA->PSubArr[proSubID]=NULL;
	SubCatB->SubLen[desSubID]=SubCatA->SubLen[proSubID];
	SubCatB->SubRank[desSubID]=SubRank;
	SubCatB->HaloChains[desSubID].HostID=HostID;
	SubCatB->HaloChains[desSubID].ProSubID=proSubID;
	memcpy(SubCatB->Property+desSubID,SubCatA->Property+proSubID,sizeof(struct SubProperty));
	//update sub_hierarchy:  need to set all to -1 for Quasis???
	if(-1==HostID)//quasi-halos,isolate them
	{
	SubCatB->sub_hierarchy[desSubID].nibs=-1;
	SubCatB->sub_hierarchy[desSubID].pre=-1;
	SubCatB->sub_hierarchy[desSubID].next=-1;
	SubCatB->sub_hierarchy[desSubID].sub=-1;
	}
	else
	{
	SubCatB->sub_hierarchy[desSubID].nibs=pro2dest[SubCatA->sub_hierarchy[proSubID].nibs];
	SubCatB->sub_hierarchy[desSubID].pre=pro2dest[SubCatA->sub_hierarchy[proSubID].pre];
	SubCatB->sub_hierarchy[desSubID].next=pro2dest[SubCatA->sub_hierarchy[proSubID].next];
	SubCatB->sub_hierarchy[desSubID].sub=pro2dest[SubCatA->sub_hierarchy[proSubID].sub];
	}
}
void migrate_src(SRCCATALOGUE *SrcCatA, SRCCATALOGUE *SrcCatB,HBTInt proSubID,HBTInt desSubID)
{
	SrcCatB->SubLen[desSubID]=SrcCatA->SubLen[proSubID];
	SrcCatB->SubLen2[desSubID]=SrcCatA->SubLen2[proSubID];
	SrcCatB->CoreFrac[desSubID]=SrcCatA->CoreFrac[proSubID];
	SrcCatB->PSubArr[desSubID]=SrcCatA->PSubArr[proSubID];SrcCatA->PSubArr[proSubID]=NULL;
	SrcCatB->PSubArr2[desSubID]=SrcCatA->PSubArr2[proSubID];SrcCatA->PSubArr2[proSubID]=NULL;
}
void mask_mainsub(CATALOGUE *CatB,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt desID,HBTInt son)
{/* CatB.HaloMask must have been initialized with all ones before calling
	* SubCatB's nonmain subs must have been filled before doing the mask
	* */
	HBTInt i,j,pid,subid,desSubID;
	HBTInt *GrpPIDs,*MainPIDs;
	desSubID=SubCatB->GrpOffset_Sub[desID];
	j=0;
	for(i=1;i<SubCatB->GrpLen_Sub[desID];i++)
	{
		subid=desSubID+i;
		j+=SubCatB->SubLen[subid];
		for(pid=0;pid<SubCatB->SubLen[subid];pid++)
			CatB->HaloMask[SubCatB->PSubArr[subid][pid]]=0;
	}
	SubCatB->SubLen[desSubID]=CatB->Len[desID]-j;
	MainPIDs=mymalloc(sizeof(HBTInt)*SubCatB->SubLen[desSubID]);
	GrpPIDs=CatB->PIDorIndex+CatB->Offset[desID];
	if(j==0)//single pro or infantry Grp,need no mask
		memcpy(MainPIDs,GrpPIDs,sizeof(HBTInt)*CatB->Len[desID]);
	else
	{
		j=0;
		for(i=0;i<CatB->Len[desID];i++)
		{
			if(CatB->HaloMask[pid=GrpPIDs[i]])
			{
				if(j==SubCatB->SubLen[desSubID])//SubLen already filled while new particles remain
				//the case when dup-particle subs occur inside this halo, so that the estimated Len is smaller than real
				MainPIDs=realloc(MainPIDs,sizeof(HBTInt)*(j+CatB->Len[desID]-i));//enlarge the memory to allow for unprocessed particles
				MainPIDs[j]=pid;
				j++;
			}
		}
		if(j>SubCatB->SubLen[desSubID])//dup-particle case, memory has been enlarged
		{
		SubCatB->SubLen[desSubID]=j;
		MainPIDs=realloc(MainPIDs,sizeof(HBTInt)*j);//adjust the memory size to fit
		}
		else if(j<SubCatB->SubLen[desSubID])//particles masked-out by other halos? strictly not expected since split_srccat has restricted all subs's partilces to be inside their host halo
		{
			fprintf(logfile,"error: Mask Fof len mismatch! \n for desID="HBTIFMT"\n,remained="HBTIFMT",expected="HBTIFMT"",desID,j,SubCatB->SubLen[desSubID]);fflush(logfile);
			exit(1);	
		}	
	}
	SubCatB->PSubArr[desSubID]=MainPIDs;
	SubCatB->SubRank[desSubID]=0;
	SubCatB->HaloChains[desSubID].ProSubID=proSubID;
	SubCatB->HaloChains[desSubID].HostID=desID;
	SubCatB->sub_hierarchy[desSubID].nibs=-1;
	SubCatB->sub_hierarchy[desSubID].pre=-1;
	SubCatB->sub_hierarchy[desSubID].next=-1;
	SubCatB->sub_hierarchy[desSubID].sub=son;
}
void mask_mainsrc(CATALOGUE *CatB,SRCCATALOGUE *SrcCatB,HBTInt desID,HBTInt desSubID,HBTInt GrpLen_Sub)
{/* CatB.HaloMaskSrc must have been initialized with all ones before calling
 */
	HBTInt i,j,pid,subid;
	HBTInt *GrpPIDs,*MainPIDs;
	MainPIDs=mymalloc(sizeof(HBTInt)*CatB->Len[desID]);
	GrpPIDs=CatB->PIDorIndex+CatB->Offset[desID];
	if(GrpLen_Sub<2)//infantry or single pro
	{
		memcpy(MainPIDs,GrpPIDs,sizeof(HBTInt)*CatB->Len[desID]);
		SrcCatB->SubLen[desSubID]=CatB->Len[desID];
		SrcCatB->PSubArr[desSubID]=MainPIDs;
	}
	else
	{
		for(i=1;i<GrpLen_Sub;i++)
		{
			subid=desSubID+i;
			for(pid=0;pid<SrcCatB->SubLen[subid];pid++)
				CatB->HaloMaskSrc[SrcCatB->PSubArr[subid][pid]]=0;
		}
		j=0;
		for(i=0;i<CatB->Len[desID];i++)
		{
			if(CatB->HaloMaskSrc[pid=GrpPIDs[i]])
			{
				MainPIDs[j]=pid;
				j++;
			}
		}
		SrcCatB->SubLen[desSubID]=j;
		SrcCatB->PSubArr[desSubID]=realloc(MainPIDs,sizeof(HBTInt)*j);
	}
	SrcCatB->SubLen2[desSubID]=-1;
	SrcCatB->CoreFrac[desSubID]=CoreFrac0;
}
void restore_mainsub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt desSubID)
{
	SubCatB->PSubArr[desSubID]=SubCatA->PSubArr[proSubID];SubCatA->PSubArr[proSubID]=NULL;
	SubCatB->SubLen[desSubID]=SubCatA->SubLen[proSubID];
	memcpy(SubCatB->Property+desSubID,SubCatA->Property+proSubID,sizeof(struct SubProperty));
}

void narrow_srccat(SRCCATALOGUE *SrcCat,SUBCATALOGUE *SubCat, HBTInt subid)
{/*those without Grp2 will have SubLen2=-1 and CoreFrac=CoreFrac0 and PSubArr2=NULL*/
	if(SubCat->SubLen[subid]==0)//unbind to death
	{
		SrcCat->SubLen[subid]=0;
		SrcCat->SubLen2[subid]=0;
		if(SrcCat->PSubArr[subid]!=NULL){free(SrcCat->PSubArr[subid]);SrcCat->PSubArr[subid]=NULL;}
		if(SrcCat->PSubArr2[subid]!=NULL){free(SrcCat->PSubArr2[subid]);SrcCat->PSubArr2[subid]=NULL;}
		SrcCat->CoreFrac[subid]=CoreFrac0;
	}
	else if(SubCat->SubLen[subid]>SrcCat->SubLen[subid])//over-growth or merger, update Src and reset Src2
	{
		if(SubCat->SubLen[subid]>1.5*SrcCat->SubLen[subid]&&SubCat->SubLen[subid]>400)//record major merger (m/M>1/2)
		fprintf(logfile,""HBTIFMT":("HBTIFMT", "HBTIFMT"),"HBTIFMT"\n",subid,SrcCat->SubLen[subid],SrcCat->SubLen2[subid],SubCat->SubLen[subid]);
		SrcCat->SubLen[subid]=SubCat->SubLen[subid];
		free(SrcCat->PSubArr[subid]);
		SrcCat->PSubArr[subid]=mymalloc(sizeof(HBTInt)*SubCat->SubLen[subid]);
		memcpy(SrcCat->PSubArr[subid],SubCat->PSubArr[subid],sizeof(HBTInt)*SubCat->SubLen[subid]);
		SrcCat->SubLen2[subid]=-1;
		myfree(SrcCat->PSubArr2[subid]);
		SrcCat->PSubArr2[subid]=NULL;
		SrcCat->CoreFrac[subid]=CoreFrac0;
	}	
	//smooth accretion or stripping
	else if(SrcCat->SubLen2[subid]>=0)  //already have Grp2
	{
		if(SubCat->SubLen[subid]>SrcCat->SubLen2[subid]) //update GrpPIDs2 when subhalo has grown big enough
		{
			SrcCat->SubLen2[subid]=SubCat->SubLen[subid];
			free(SrcCat->PSubArr2[subid]);
			SrcCat->PSubArr2[subid]=mymalloc(sizeof(HBTInt)*SrcCat->SubLen2[subid]);
			memcpy(SrcCat->PSubArr2[subid],SubCat->PSubArr[subid],sizeof(HBTInt)*SubCat->SubLen[subid]);
			SrcCat->CoreFrac[subid]=((HBTReal)SrcCat->SubLen2[subid])/((HBTReal)SrcCat->SubLen[subid])/MassRelax_Input;//?? necessary to change this?
		}	
		else if(SubCat->SubLen[subid]<SrcCat->SubLen2[subid]/MassRelax_Input)//update GrpPIDs with GrpPIDs2, and update GrpPIDs2 with the current sub (or its nearest progenitor?)
		{
			SrcCat->SubLen[subid]=SrcCat->SubLen2[subid];
			SrcCat->SubLen2[subid]=SubCat->SubLen[subid];
			free(SrcCat->PSubArr[subid]);
			SrcCat->PSubArr[subid]=SrcCat->PSubArr2[subid];
			SrcCat->PSubArr2[subid]=mymalloc(sizeof(HBTInt)*SrcCat->SubLen2[subid]);			
			memcpy(SrcCat->PSubArr2[subid],SubCat->PSubArr[subid],sizeof(HBTInt)*SubCat->SubLen[subid]);
			SrcCat->CoreFrac[subid]=((HBTReal)SrcCat->SubLen2[subid])/((HBTReal)SrcCat->SubLen[subid])/MassRelax_Input;/*every time when GrpLen2 changes, 
																											* update CoreFrac to ensure we will be 
																											* operating in the same range ,
																											* that is, the lower limit sub will be unbind() with
																											* near real CoreFrac, rather than a possibly much higher
																											* CoreFrac if not updated*/
		}
	}
	//no Grp2 yet
	else if(SubCat->SubLen[subid]<=SrcCat->SubLen[subid]/MassRelax_Input)// register GrpPIDs2 when subhalo has reduced enough, and adapt CoreFrac
	{
		SrcCat->SubLen2[subid]=SubCat->SubLen[subid];
		SrcCat->PSubArr2[subid]=mymalloc(sizeof(HBTInt)*SubCat->SubLen[subid]);//need not free first, GrpPIDs2 haven't been allocated when Len2==Len;
		memcpy(SrcCat->PSubArr2[subid],SubCat->PSubArr[subid],sizeof(HBTInt)*SubCat->SubLen[subid]);
		SrcCat->CoreFrac[subid]=((HBTReal)SrcCat->SubLen2[subid])/((HBTReal)SrcCat->SubLen[subid])/MassRelax_Input;//?? necessary to change this?
	}
}

void PARAinit_mask(CATALOGUE *Cat,HBTInt mtype)
{
	HBTInt i;
	if(mtype==0)//init_both
	{
		#pragma omp single
		{
		Cat->HaloMask=mymalloc(sizeof(char)*NP_DM);
		Cat->HaloMaskSrc=mymalloc(sizeof(char)*NP_DM);
		}
		#pragma omp for 
		for(i=0;i<NP_DM;i++)
		{
			Cat->HaloMask[i]=1;
			Cat->HaloMaskSrc[i]=1;
		}
	}
	else if(mtype==1)			//init sub
	{
		#pragma omp single
		Cat->HaloMask=mymalloc(sizeof(char)*NP_DM);
		#pragma omp for
		for(i=0;i<NP_DM;i++)
			Cat->HaloMask[i]=1;
	}
	else if(mtype==2)		//init src
	{
		#pragma omp single
		Cat->HaloMaskSrc=mymalloc(sizeof(char)*NP_DM);
		#pragma omp for
		for(i=0;i<NP_DM;i++)
			Cat->HaloMaskSrc[i]=1;
	}
	else
	#pragma omp single
	{
		fprintf(logfile,"error using init_mask, wrong mtype\n");fflush(logfile);
		exit(1);
	}
}

void mask_src_recursive(HBTInt subid,SUBCATALOGUE *SubCat, SRCCATALOGUE *SrcCat,char *HaloMaskSrc)
{
	HBTInt i,pid,son,sib;
	HBTInt *SubArr, Len;
	
	if((son=SubCat->sub_hierarchy[subid].sub)>=0)//has sons,mask them recursively to update its reservior
	{	
		sib=son;
		while(sib>=0)
		{
			mask_src_recursive(sib,SubCat,SrcCat,HaloMaskSrc);
			sib=SubCat->sub_hierarchy[sib].next;
		}
	}
	//mask this sub
	if(SrcCat->SubLen[subid])
	{
	SubArr=mymalloc(sizeof(HBTInt)*SrcCat->SubLen[subid]);
	Len=0;
	for(i=0;i<SrcCat->SubLen[subid];i++)
	{
		pid=SrcCat->PSubArr[subid][i];
		if(HaloMaskSrc[pid])//still not masked-out by its sons, take this
		{
			SubArr[Len]=pid;
			Len++;
			HaloMaskSrc[pid]=0;//mask it out
		}
	}
	free(SrcCat->PSubArr[subid]);
	SrcCat->PSubArr[subid]=realloc(SubArr,sizeof(HBTInt)*Len);
	SrcCat->SubLen[subid]=Len;
	}
}

static int compare_host(const void *a, const void *b)//used to sort desendent id in ascending order
{
  if(((struct Chain_data *) a)->HostID > ((struct Chain_data *) b)->HostID)
    return +1;

  if(((struct Chain_data *) a)->HostID < ((struct Chain_data *) b)->HostID)
    return -1;

  return 0;
}
static SUBCATALOGUE *SubCatCMP;//this is used as a bridge var between init_dessub()  and compare_proLen() for sorting
static int compare_proLen(const void *a, const void *b)//used to sort progenitor size in descending order
	{
		  //~ extern SUBCATALOGUE *SubCatCMP;//some caution here.........
		  if(SubCatCMP->SubLen[((struct Chain_data *) a)->ProSubID] > SubCatCMP->SubLen[((struct Chain_data *) b)->ProSubID] )
			return -1;

		 if(SubCatCMP->SubLen[((struct Chain_data *) a)->ProSubID] < SubCatCMP->SubLen[((struct Chain_data *) b)->ProSubID] )
			return +1;

		  return 0;
	}
void PARAinit_dessub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB, struct LinkInfo *linkinfo)
{/*output: SubCatA, SubCatB, linkinfo *
	*purpose:
	* 				compute statistics: birth,death,quasi,and linkinfo (procount,chainoffset,pro2dest)
	* 				initialize SubCatB with the statistics
	* 				update SubCatA.sub_hierarchy and sort SubCatA.HaloChains
	* 		*/
	
	HBTInt i,desID,subhaloid,subhaloid2,son;
	HBTInt BirthCount,DeathCount,GrpLen_Sub,GrpOffset_Sub,GrpOffset_chain;
	struct Chain_data *HaloChains_alive,*GrpChain;
	
	#pragma omp single copyprivate(HaloChains_alive)
	{
	BirthCount=0;
	DeathCount=0;
	GrpLen_Sub=0;
	GrpOffset_Sub=0;
	GrpOffset_chain=0;
	
	SubCatCMP=SubCatA;
		
	linkinfo->ProCount=calloc(SubCatB->Ngroups+1,sizeof(HBTInt));linkinfo->ProCount++;
	linkinfo->ChainOffset=calloc(SubCatB->Ngroups+1,sizeof(HBTInt));linkinfo->ChainOffset++;
	for(i=0;i<SubCatA->Nsubs;i++)
	{
		if(SubCatA->SubLen[i])
			linkinfo->ProCount[SubCatA->HaloChains[i].HostID]++;
		else
		{
			SubCatA->HaloChains[i].HostID=-2;//death
			DeathCount++;
		}
	}
	qsort(SubCatA->HaloChains,SubCatA->Nsubs,sizeof(struct Chain_data),compare_host);//consider improving efficiency by manual sorting here
	SubCatB->GrpOffset_Sub=mymalloc(sizeof(HBTInt)*SubCatB->Ngroups);   SubCatB->GrpLen_Sub=mymalloc(sizeof(HBTInt)*SubCatB->Ngroups);
	/*=====prepare Len and Offset info for parallel access to halo; get birthcount to decide Nsubs===========*/
	for(desID=0;desID<SubCatB->Ngroups;desID++)//we don't count background(desID=-1) as birth in any case
	{		
		linkinfo->ChainOffset[desID]=GrpOffset_chain;
		if((GrpLen_Sub=linkinfo->ProCount[desID]))
		{			
			GrpOffset_chain+=GrpLen_Sub;
		}
		else //no progenitor
		{
			GrpLen_Sub=1;
			BirthCount++;//====what if this baby dies immediately?Then this sub would be of 0 length but still ocuppies one position
		}
		SubCatB->GrpLen_Sub[desID]=GrpLen_Sub;
		SubCatB->GrpOffset_Sub[desID]=GrpOffset_Sub;
		GrpOffset_Sub+=GrpLen_Sub;
	}
	linkinfo->ChainOffset[-1]=-linkinfo->ProCount[-1];//the background offset,i.e,quasi-halos' offset
	
	SubCatB->Nbirth=BirthCount;    SubCatB->Ndeath=DeathCount;       SubCatB->NQuasi=linkinfo->ProCount[-1]; SubCatB->Nsplitter=SubCatA->Nsplitter;
	SubCatB->Nsubs=SubCatA->Nsubs-DeathCount+BirthCount;
	create_sub_cat(SubCatB);
	fprintf(logfile,"Nsublast="HBTIFMT"\tNsub="HBTIFMT"\tBirth="HBTIFMT"\tDeath="HBTIFMT"\tQuasi="HBTIFMT"\tSplitter="HBTIFMT"\n",SubCatA->Nsubs,SubCatB->Nsubs,BirthCount,DeathCount,SubCatB->NQuasi,SubCatB->Nsplitter);//here Nsublast includes Nsplitters
	fflush(logfile);

	/*construct pro2dest table for updates of sub_hierarchy; also find mainsub and update SubCatTmp.sub_hierarchy with old subid
		*	be aware sub_hierarchy have to be cleaned first, i.e, sub_cross marked and dead sub kicked, 
		*	so that after this step we only need to replace old subid with dessubid to finish sub_hierarchy */
	HaloChains_alive=SubCatA->HaloChains+DeathCount+SubCatB->NQuasi;
	linkinfo->pro2dest=mymalloc(sizeof(HBTInt)*(SubCatA->Nsubs+1));//remember to free it; Nsubs+1 to account for id=-1 which means nobody
	linkinfo->pro2dest[0]=-1; linkinfo->pro2dest+=1;
	for(i=0;i<DeathCount;i++)
		linkinfo->pro2dest[SubCatA->HaloChains[i].ProSubID]=-1;
	GrpChain=SubCatA->HaloChains+DeathCount;//skip death subs, point to quasi now
	GrpLen_Sub=SubCatB->NQuasi;
	if(GrpLen_Sub>1)
		qsort(GrpChain,GrpLen_Sub,sizeof(struct Chain_data),compare_proLen);
	for(i=GrpLen_Sub-1;i>=0;i--)//store (and hence lable) quasi-halos from least massive to most massive, different from normal groups
		linkinfo->pro2dest[GrpChain[i].ProSubID]=SubCatB->Nsubs-i-1;//consider modify this to normally sort quasi-halos
	}
	#pragma omp for schedule(dynamic,1)
	for(desID=0;desID<SubCatB->Ngroups;desID++)
	{
		if((GrpLen_Sub=linkinfo->ProCount[desID]))
		{
			GrpChain=HaloChains_alive+linkinfo->ChainOffset[desID];
			if(GrpLen_Sub>1)
			{
				qsort(GrpChain,GrpLen_Sub,sizeof(struct Chain_data),compare_proLen);
				break_out_sub_family(GrpChain[0].ProSubID,SubCatA);//mark mainsub as crossing-out
			}
			subhaloid=GrpChain[0].ProSubID;
			linkinfo->pro2dest[subhaloid]=SubCatB->GrpOffset_Sub[desID];
			son=SubCatA->sub_hierarchy[subhaloid].sub;
			while(son>=0)
			{
				SubCatA->sub_hierarchy[subhaloid].next=son;
				son=SubCatA->sub_hierarchy[son].next;
			}//store the last son temporarily	
			for(i=1;i<GrpLen_Sub;i++)
			{
				subhaloid2=GrpChain[i].ProSubID;
				linkinfo->pro2dest[subhaloid2]=SubCatB->GrpOffset_Sub[desID]+i;
				if(SubCatA->sub_hierarchy[subhaloid2].nibs<0)//a highest-level sub
				{
					if((son=SubCatA->sub_hierarchy[subhaloid].next)<0)//mainsub has no son yet
					{
						SubCatA->sub_hierarchy[subhaloid].sub=subhaloid2;
						SubCatA->sub_hierarchy[subhaloid].next=subhaloid2;
					}
					else//already has sons
					{
						SubCatA->sub_hierarchy[son].next=subhaloid2;
						SubCatA->sub_hierarchy[subhaloid2].pre=son;
						SubCatA->sub_hierarchy[subhaloid].next=subhaloid2;
					}
					SubCatA->sub_hierarchy[subhaloid2].nibs=subhaloid;
				}
			}
			SubCatA->sub_hierarchy[subhaloid].next=-1;
		}
	}
}

/*===========below are debugging functions============*/
static void errorfun(HBTInt type,HBTInt desID,HBTInt i,HBTInt pid)
{
	printf("type"HBTIFMT", grp"HBTIFMT",p"HBTIFMT",pid"HBTIFMT"\n",type,desID,i,pid);
}
void mask_mainsub_check(CATALOGUE *CatB,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt desID,HBTInt son)
{/* CatB.HaloMask must have been initialized with all ones before calling
	* SubCatB's nonmain subs must have been filled before doing the mask
	* */
	HBTInt i,j,k,pid,subid,desSubID;
	HBTInt *GrpPIDs,*MainPIDs;
	short *catmask;
	catmask=calloc(NP_DM,sizeof(short));
	desSubID=SubCatB->GrpOffset_Sub[desID];
	GrpPIDs=CatB->PIDorIndex+CatB->Offset[desID];
	for(i=0;i<CatB->Len[desID];i++)
	{
		if(CatB->HaloMask[GrpPIDs[i]]==0)
			errorfun(0,desID,i,GrpPIDs[i]);
		catmask[GrpPIDs[i]]=1;
	}
	j=0;
	for(i=1;i<SubCatB->GrpLen_Sub[desID];i++)
	{
		subid=desSubID+i;
		j+=SubCatB->SubLen[subid];
		for(pid=0;pid<SubCatB->SubLen[subid];pid++)
		{
			if(CatB->HaloMask[SubCatB->PSubArr[subid][pid]]==0)
			errorfun(1,subid,i,SubCatB->PSubArr[subid][pid]);
			else 
			CatB->HaloMask[SubCatB->PSubArr[subid][pid]]=0;
			if(catmask[SubCatB->PSubArr[subid][pid]]==0)//outlier
			printf("subid"HBTIFMT",i"HBTIFMT",ind"HBTIFMT",pid"HBTIFMT"\n",subid,i,pid,SubCatB->PSubArr[subid][pid]);
		}
	}
	SubCatB->SubLen[desSubID]=CatB->Len[desID]-j;
	MainPIDs=mymalloc(sizeof(HBTInt)*SubCatB->SubLen[desSubID]);
	
	if(j==0)//single pro or infantry Grp,need no mask
		memcpy(MainPIDs,GrpPIDs,sizeof(HBTInt)*CatB->Len[desID]);
	else
	{
		j=0;
		for(i=0;i<CatB->Len[desID];i++)
		{
			if(CatB->HaloMask[pid=GrpPIDs[i]])
			{
				MainPIDs[j]=pid;
				j++;
			}
		}
		if(j!=SubCatB->SubLen[desSubID])//what about if there're duplicate particles???
		{
			fprintf(logfile,"error: Mask Fof len mismatch! \n for desID="HBTIFMT"\n,remained="HBTIFMT",expected="HBTIFMT"",desID,j,SubCatB->SubLen[desSubID]);fflush(logfile);
			exit(1);
		}
	}
	SubCatB->PSubArr[desSubID]=MainPIDs;
	SubCatB->SubRank[desSubID]=0;
	SubCatB->HaloChains[desSubID].ProSubID=proSubID;
	SubCatB->HaloChains[desSubID].HostID=desID;
	SubCatB->sub_hierarchy[desSubID].nibs=-1;
	SubCatB->sub_hierarchy[desSubID].pre=-1;
	SubCatB->sub_hierarchy[desSubID].next=-1;
	SubCatB->sub_hierarchy[desSubID].sub=son;
}
