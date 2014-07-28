#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
char outputdir[1024]=SUBCAT_DIR;
char fofdir[1024]=GRPCAT_DIR;
char snapdir[1024]=SNAPSHOT_DIR;
	
HBTInt SnapshotNum=0,SnapRange[2]={0};
CATALOGUE CatB;
SUBCATALOGUE SubCatA,SubCatB;
SRCCATALOGUE SrcCatA,SrcCatB;

struct LinkInfo linkinfo;
struct Chain_data * HaloChains_alive, *GrpChain;//HaloChains_alive is offseted by DeathCount and Nquasi to contain only relevant chains for CatB
HBTInt Lmain_removed,*PInd_main_removed;
HBTInt i,desID,haloid,subhaloid,GrpLen_Sub;
HBTInt Nff,Nfm;//NfailedFoF: failed fof withou progenitor ,  NfailedMain: failed mainsub with one or more pro's
char buf[1024],filemode[4];
time_t time_start,time_end,tward[20]; //for program timing

FILE *fptime;

time_start=time(NULL);
#ifdef _OPEN_MP
#ifdef HALO_PARA
omp_set_nested(0);
#endif
#endif
/*=================program initialization==================*/
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
	 
	 if(SnapRange[0]>IniSnap)//add logs to existing file
			sprintf(filemode,"%s","a");
	else//create new logfiles
			sprintf(filemode,"%s","w");
	
	mkdir(SUBCAT_DIR,0755);		
	if(0==(strcmp(LOGFILE_NAME,"stdout"))) 
		logfile=stdout;
	else
	{	
		sprintf(buf,"%s/%s",outputdir,LOGFILE_NAME);
		if((logfile=fopen(buf,filemode))==NULL)	{fprintf(stderr,"Error opening file %s\n",buf);exit(1);}
	}
	
	sprintf(buf,"%s/timing.log",SUBCAT_DIR);
	myfopen(fptime,buf,filemode);

	sprintf(buf,"%s/splitters",outputdir); //sp2pro table
	mkdir(buf,0755);
	sprintf(buf,"%s/pro2dest",outputdir); //pro2dest table
	mkdir(buf,0755);
	
	FILE * VERFILE;
	sprintf(buf,"%s/VER%s",outputdir,getenv("HBT_VERSION"));
	if((VERFILE=fopen(buf,"a"))==NULL)	{fprintf(stderr,"Error opening file %s\n",buf);exit(1);}
	fclose(VERFILE);
#ifndef UNBIND_FOF		
	if(SnapRange[0]>IniSnap)//pick up a previous subcat to continue,assign snapnum first
	{		
			fprintf(logfile,"\nrestarting program %s at %s from Snapshot "HBTIFMT"\n",argv[0],ctime(&time_start),SnapRange[0]);
			SnapshotNum=SnapRange[0]-1;
			load_sub_catalogue(SnapshotNum,&SubCatB,outputdir);
			load_src_catalogue(SnapshotNum,&SrcCatB,outputdir);
	}
	else//forge SnapshotNum=-1
#endif	  
	{	
			SubCatB.Ngroups=0;SubCatB.GrpOffset_Sub=NULL;SubCatB.GrpLen_Sub=NULL;
			SubCatB.Nsubs=0;	create_sub_cat(&SubCatB);
			SrcCatB.Nsubs=0;   create_src_cat(&SrcCatB);
	}	
			SubCatA.Ngroups=0;SubCatA.GrpOffset_Sub=NULL;SubCatA.GrpLen_Sub=NULL;
			SubCatA.Nsubs=0;   create_sub_cat(&SubCatA);
			SrcCatA.Nsubs=0;    create_src_cat(&SrcCatA);
}
/*=======================end initialization==========================*/

for(SnapshotNum=SnapRange[0];SnapshotNum<=SnapRange[1];SnapshotNum++)
{
	tward[0]=time(NULL);
	load_particle_data(SnapshotNum,snapdir);
	load_group_catalogue(SnapshotNum,&CatB,fofdir);	
	tward[1]=time(NULL);
	//now convert PID to particle address
	fill_PIDHash();
	fresh_ID2Index(&CatB,FRSH_GRPCAT); 
	fresh_ID2Index(&SrcCatB,FRSH_SRCCAT);//now CatB.PIDorIndex has been occupied with Index to the freshly loaded Pdat
								//also refresh subPindex according to the new Pdat.
	free_PIDHash();
	tward[2]=time(NULL);
	prepare_ind2halo(&CatB);															
	
#ifdef HALO_PARA	
  #pragma omp parallel
#endif
  {//\\begin omp parallel	
	PARAsplit_srccat(&CatB,&SrcCatB,SnapshotNum);
	#ifdef HALO_PARA
	#pragma omp master
	#endif
	tward[3]=time(NULL);
	#ifdef HALO_PARA
	#pragma omp single
	#endif
	fprintf(logfile,"DeathSrc="HBTIFMT"\n",SrcCatB.NDeathSp);fflush(logfile);
	
	PARAmake_srcsub(&SubCatA,&SubCatB,&SrcCatA,&SrcCatB);
	/*=update sub_hierarchy to mark cross-out and kick out dead sub , 
	 * this is necessary since we have splitted srccat and finished desHalo finding
	 *Note this just break the new quasi-halos from their old hosts,
	 * but not fully break the sub_hierarchy among those quasi-halos which come out together to the background.
	 * since old quasi-halos are all isolated,no need to markcross for them =*/
	#ifdef HALO_PARA
	#pragma omp for schedule(dynamic,1)
	#endif
	for(haloid=0;haloid<SubCatA.Ngroups;haloid++)
	{
		markcross_N_kickdeath_recursive(SubCatA.GrpOffset_Sub[haloid],&SubCatA);
	}
	#ifdef HALO_PARA
	#pragma omp master
	#endif
	tward[4]=time(NULL);
//till now the srccat will be unbinded in its own des halo, so the dessubs will have no particle from other halos,safe for mask_mainsub later
	//unbind the last subcat in the current epoch
	#ifdef HALO_PARA
	#pragma omp for private(subhaloid,Lmain_removed,PInd_main_removed) schedule(dynamic,1)
	#endif
	#ifdef SAT_ACCR_ON
	for(subhaloid=0;subhaloid<SubCatA.Nsubs;subhaloid++)
	{
		if(SubCatA.sub_hierarchy[subhaloid].nibs<0&&SubCatA.SubLen[subhaloid]>0)//unbind those living head-subs
		{
		unbind_sub_recursive(subhaloid,&Lmain_removed,&PInd_main_removed,&SubCatA,&SrcCatA);
		myfree(PInd_main_removed);//PInd_main_removed can be NULL if no src or all bound
		}
	}
	#else
		for(subhaloid=0;subhaloid<SubCatA.Nsubs;subhaloid++)
	{
		if(SubCatA.SubLen[subhaloid]>0)//unbind those living subs
		{
		unbind(SubCatA.SubLen+subhaloid,SubCatA.PSubArr+subhaloid,SubCatA.Property+subhaloid,&Lmain_removed,&PInd_main_removed,SrcCatA.CoreFrac[subhaloid]);
		narrow_srccat(&SrcCatA,&SubCatA, subhaloid);
		myfree(PInd_main_removed);//PInd_main_removed can be NULL if no src or all bound
		}
	}
	#endif
	#ifdef HALO_PARA
	#pragma omp master
	#endif
	tward[5]=time(NULL);
	// kick out dead sub due to unbinding
	#ifdef HALO_PARA
	#pragma omp for schedule(dynamic,1)
	#endif
	for(subhaloid=0;subhaloid<SubCatA.Nsubs;subhaloid++)
	{		
		if(SubCatA.sub_hierarchy[subhaloid].nibs<0)//take head-subs
			kickdeath_recursive(subhaloid,&SubCatA);
	}
	/*mask srccat to eliminate duplicate particles; 
	*since we use relaxed srccat and individual update of srcsub, plus recursive unbinding
	*  it is possible to have duplicate particles, 
	* (also possible is particles from other halos, so it is not guaranteed that mask_mainsrc will be inside one halo) 
	* and masking_src here ensures the srccat have no dup particles later (parallelization may permit, in quite rare case of data-racing, dup-particles though).
	* however, mask_mainsrc later would restore the masked-out particles in crossed-out or splitted-out single-sub halos, so duplicate particles remains,
	* and those dup particles only get to be eliminated when the dup-particle subs merge-in again,at which time the sub catalogue would contain dup-particle sub inside one halo
	*since this is done along subhaloid sequence, the bigger halo in prev snap will have priority to retain the dup-particles
	 * and the case is similar in mask_mainsrc() later*/
	/*improvements needed to make the dup-particles claimed by the sub in those particles' host-halo,e.g., by the splitted or crossed-out part*/
	PARAinit_mask(&CatB,2);//init_src_mask
	#ifdef HALO_PARA
	#pragma omp for schedule(dynamic,1) //when parallelized, dup-particles can in general be retained in bigger halos but not strictly the case
	#endif
	for(subhaloid=0;subhaloid<SubCatA.Nsubs;subhaloid++)
	{
		if(SubCatA.sub_hierarchy[subhaloid].nibs<0)//take head-subs
			mask_src_recursive(subhaloid,&SubCatA,&SrcCatA,CatB.HaloMaskSrc);
	}
	#ifdef HALO_PARA
	#pragma omp single		
	#endif	
	{
	free(CatB.HaloMaskSrc);
	SubCatB.Ngroups=CatB.Ngroups;	
	}
	PARAinit_dessub(&SubCatA,&SubCatB,&linkinfo);
	PARAinit_mask(&CatB,0);	
	#ifdef HALO_PARA
	#pragma omp single		
	#endif	
	{
	SrcCatB.Nsubs=SubCatB.Nsubs;	create_src_cat(&SrcCatB);
	SrcCatB.NDeathSp=SrcCatA.NDeathSp;
	Nff=0;Nfm=0;
	HaloChains_alive=SubCatA.HaloChains+SubCatB.Ndeath+SubCatB.NQuasi;
    }
	#ifdef HALO_PARA
	#pragma omp master
	#endif
	tward[6]=time(NULL);
	#ifdef HALO_PARA
	#pragma omp for private(desID,GrpChain,subhaloid,GrpLen_Sub,i,Lmain_removed,PInd_main_removed) reduction(+:Nff,Nfm) schedule(dynamic,1)
	#endif
	for(desID=0;desID<CatB.Ngroups;desID++)//process FoF_dest
	{
		GrpChain=HaloChains_alive+linkinfo.ChainOffset[desID];
		subhaloid=SubCatB.GrpOffset_Sub[desID];
		if((GrpLen_Sub=linkinfo.ProCount[desID]))//ok,need to mask to find main-sub
		{
			/*========transfer non-main subs and mask particles to update main sub======*/
			for(i=1;i<GrpLen_Sub;i++)
			{
				migrate_sub(&SubCatA,&SubCatB,GrpChain[i].ProSubID,i,desID,linkinfo.pro2dest);
				migrate_src(&SrcCatA,&SrcCatB,GrpChain[i].ProSubID,subhaloid+i);
			}
			mask_mainsub(&CatB,&SubCatB,GrpChain[0].ProSubID,desID,linkinfo.pro2dest[SubCatA.sub_hierarchy[GrpChain[0].ProSubID].sub]);
			if(!unbind(SubCatB.SubLen+subhaloid, SubCatB.PSubArr+subhaloid,SubCatB.Property+subhaloid,&Lmain_removed,&PInd_main_removed,CoreFrac0))
			{//so, failed to unbind the maskedfof, retreat to MainCat
				if(CatB.Len[desID]>400)
				#ifdef HALO_PARA
				#pragma omp critical (log_write)
				#endif
				fprintf(logfile,""HBTIFMT"\t",desID);
				restore_mainsub(&SubCatA,&SubCatB,GrpChain[0].ProSubID,subhaloid);
				migrate_src(&SrcCatA,&SrcCatB,GrpChain[0].ProSubID,subhaloid);
				Nfm++;				
			}
			else
			{
				mask_mainsrc(&CatB,&SrcCatB,desID,subhaloid,GrpLen_Sub);
				free(SubCatA.PSubArr[GrpChain[0].ProSubID]);
				free(SrcCatA.PSubArr[GrpChain[0].ProSubID]);
				myfree(SrcCatA.PSubArr2[GrpChain[0].ProSubID]);
			}
			myfree(PInd_main_removed);
		}
		else//infantry halo 
		{
			mask_mainsub(&CatB,&SubCatB,-1,desID,-1);
			/*============unbind the FoF=========================*/
			if(!unbind(SubCatB.SubLen+subhaloid, SubCatB.PSubArr+subhaloid,SubCatB.Property+subhaloid,&Lmain_removed,&PInd_main_removed,CoreFrac0))
				Nff++;
			myfree(PInd_main_removed);
			mask_mainsrc(&CatB,&SrcCatB,desID,subhaloid,GrpLen_Sub);//give the entire fof as src no matter unbind survives or not,since main-subs always occupy an ID
		}
	}
  }//\\end omp parallel
	tward[7]=time(NULL);
	free(CatB.HaloMask);
	free(CatB.HaloMaskSrc);
	if(0==CatB.Ngroups){ Nff=0;Nfm=0;} //in case the unexecuted parallel loop make Nff and Nfm uninitialized
	fprintf(logfile,"\nNff="HBTIFMT"\tNfm="HBTIFMT"\n",Nff,Nfm);fflush(logfile);
	/*=================sort and copy quasi-halos=======================*/
	if((GrpLen_Sub=SubCatB.NQuasi))//ok, we have quasi-halos
	{
		GrpChain=HaloChains_alive+linkinfo.ChainOffset[-1];
		#ifdef HALO_PARA
		#pragma omp parallel for schedule(dynamic,1) if(GrpLen_Sub>NParaMin) 
		#endif
		for(i=0;i<GrpLen_Sub;i++)
		{
			migrate_sub(&SubCatA,&SubCatB,GrpChain[i].ProSubID,i,-1,linkinfo.pro2dest);
			migrate_src(&SrcCatA,&SrcCatB,GrpChain[i].ProSubID,linkinfo.pro2dest[GrpChain[i].ProSubID]);
		}
	}
	complete_N_save(&SubCatB,&SrcCatB,SnapshotNum,outputdir);
	save_pro2dest(SnapshotNum-1,linkinfo.pro2dest,SubCatA.Nsubs,outputdir);
	fflush(logfile);
//	fflush(fpsp);
	
	/*=====refill_main_sub and infantry to_fof===========*/
	free_particle_data();
	free_catalogue(&CatB);
	free(linkinfo.pro2dest-1);
	free(linkinfo.ChainOffset-1);
	free(linkinfo.ProCount-1);
	tward[8]=time(NULL);
	fprintf(fptime,"%ld ",(long)SnapshotNum);
	for(i=1;i<9;i++)
	fprintf(fptime,"%ld ",tward[i]-tward[i-1]);
	fprintf(fptime,"\n");fflush(fptime);
}
fclose(fptime);
			
/*==============program cleaning-up===========*/
for(i=0;i<SubCatB.Nsubs;i++)
{
	free(SubCatB.PSubArr[i]);
	free(SrcCatB.PSubArr[i]);
	myfree(SrcCatB.PSubArr2[i]);
}
free_sub_catalogue(&SubCatB);
free_src_catalogue(&SrcCatB);
free_sub_catalogue(&SubCatA);
free_src_catalogue(&SrcCatA);

time_end=time(NULL);
fprintf(logfile,"Program %s Finished in %ld minutes (or %ld hours) (or %ld seconds).\n",argv[0], (time_end-time_start)/60,(time_end-time_start)/3600,time_end-time_start);
fclose(logfile);

//remove write permissions from all to protect the files
#ifndef UNBIND_FOF
sprintf(buf,"chmod a-w %s/* -R",SUBCAT_DIR);
system(buf);
#endif
//create additional directories for future use
sprintf(buf,"%s/profile",SUBCAT_DIR); //density profile dir, not used by HBT, to be used by post-processing
mkdir(buf,0755);
sprintf(buf,"%s/history",SUBCAT_DIR);//evolution_cat dir, to be used by post－proc
mkdir(buf,0755);
sprintf(buf,"%s/anal",SUBCAT_DIR);//analysis dir, to be used by post-proc
mkdir(buf,0755);

return 0;
}


/*==Bug Fix:
 * mask_src_recursive() with all head-subs rather than with main-sub's from prev. snapshot. 18.05.2009
 * ==*/
















