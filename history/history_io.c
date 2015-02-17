#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "history_vars.h"
#include "history_proto.h"

HBTReal load_halo_size(HALOSIZE *halosize,HBTInt Ngroups,HBTInt Nsnap)
{
	HBTInt i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,(int)Nsnap);
	myfopen(fp,buf,"r");
	IO_HEADER h;
	load_particle_header_into(Nsnap,SNAPSHOT_DIR, &h);
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,13*4L,SEEK_CUR);
		fread(&halosize[i].mass,sizeof(int),1,fp);
		fread(halosize[i].Mvir,sizeof(int),3,fp);
		fread(halosize[i].Rvir,sizeof(float),3,fp);
		fread(halosize[i].flag_badvir,sizeof(int),3,fp);
		fread(&halosize[i].flag_fakehalo,sizeof(int),1,fp);
		if(halosize[i].flag_badvir[0])
		{
			halosize[i].Rvir[0]=comoving_virial_radius_header(halosize[i].mass, &h);
			halosize[i].Mvir[0]=halosize[i].mass;
		}
	}
	fclose(fp);
	return h.time;//return current scale factor
}
HBTInt load_halo_concentration(float *halocon,HBTInt Nsnap)
{
	HBTInt grpid,Ngroups,Ngroups2,halostatus;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_param_%03d",SUBCAT_DIR,(int)Nsnap);
	myfopen(fp,buf,"r");
	fread(&Ngroups,sizeof(int),1,fp);
	for(grpid=0;grpid<Ngroups;grpid++)
	{
		fread(&halostatus,sizeof(int),1,fp);
		fseek(fp,4*sizeof(float),SEEK_CUR);
		//~ fread(halostatus+grpid,sizeof(int),1,fp);
		fread(halocon+grpid,sizeof(float),1,fp);
		if(halostatus!=0) halocon[grpid]=-1; //set concentration to -1 if fitting not successful
		fseek(fp,sizeof(float)*2,SEEK_CUR);
	}	
	fread(&Ngroups2,sizeof(int),1,fp);
	if(Ngroups2!=Ngroups)
	{
		printf("error:Ngroups="HBTIFMT","HBTIFMT" do not match when loading \n %s\n" 
			"probably file corruption or different file format\n",
			Ngroups,Ngroups2,buf);
		exit(-1);
	}
	fclose(fp);
	return Ngroups;
}
HBTInt read_mainsubid(HBTInt Nsnap,HBTInt hostid)
{
FILE *fd;
char buf[1024];
HBTInt subid,Ngroups;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(HBTInt), 1, fd);
  if(hostid<0||hostid>=Ngroups)
  {
	  printf("error: wrong hostid to read for mainsub\n");
	  exit(1);
  }
  fseek(fd,sizeof(HBTInt)*(2+Ngroups+hostid),SEEK_CUR);
  fread(&subid,sizeof(HBTInt),1,fd);
  fclose(fd);
  return subid;
}

HBTInt read_Ngroups(HBTInt Nsnap)
{
FILE *fd;
char buf[1024];
HBTInt Ngroups;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(HBTInt), 1, fd);
  
  fclose(fd);
  return Ngroups;
}
HBTInt load_mainsubid(HBTInt Nsnap, HBTInt *(*MainSubID))
{
  FILE *fd;
char buf[1024];
HBTInt subid,Ngroups,Nsubs,Nids;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Nsubs,sizeof(HBTInt),1,fd);
  fread(&Nids, sizeof(HBTInt), 1, fd);
  fseek(fd,sizeof(HBTInt)*Ngroups,SEEK_CUR);//skip grplen_sub
  *MainSubID=mymalloc(sizeof(HBTInt)*Ngroups);
  fread(*MainSubID, sizeof(HBTInt), Ngroups, fd);//grpoffset_sub
  fclose(fd);
  return Ngroups;
}
HBTInt read_subpos(HBTInt Nsnap,HBTxyz *(*snappos))
{
FILE *fd;
char buf[1024];
HBTInt subid,Ngroups,Nsubs;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Nsubs,sizeof(HBTInt),1,fd);
  *snappos=mymalloc(sizeof(HBTReal)*3*Nsubs);
  fseek(fd,sizeof(HBTInt)*(1+Ngroups*2+Nsubs*3),SEEK_CUR);
  fseek(fd,(sizeof(struct Chain_data)+sizeof(struct Hierarchy))*Nsubs,SEEK_CUR);
  for(subid=0;subid<Nsubs;subid++)
  {
  fread(*snappos+subid,sizeof(HBTReal),3,fd);
  fseek(fd,sizeof(HBTReal)*8,SEEK_CUR);
  }
  fclose(fd);
  return Nsubs;
}
HBTInt read_subvel(HBTInt Nsnap,HBTxyz *(*snapvel))
{
FILE *fd;
char buf[1024];
HBTInt subid,Ngroups,Nsubs;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Nsubs,sizeof(HBTInt),1,fd);
  *snapvel=mymalloc(sizeof(HBTReal)*3*Nsubs);
  fseek(fd,sizeof(HBTInt)*(1+Ngroups*2+Nsubs*3),SEEK_CUR);
  fseek(fd,(sizeof(struct Chain_data)+sizeof(struct Hierarchy))*Nsubs,SEEK_CUR);
  fseek(fd,sizeof(HBTReal)*3,SEEK_CUR);
  for(subid=0;subid<Nsubs;subid++)
  {
  fread(*snapvel+subid,sizeof(HBTReal),3,fd);
  fseek(fd,sizeof(HBTReal)*8,SEEK_CUR);
  }
  fclose(fd);
  return Nsubs;
}
size_t load_brfcat(HBTInt Nsnap,BRFCAT *Cat)
{
	SUBCATALOGUE SubCat;
	HBTInt *pro2dest,Npro,i;
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	Cat->Nsnap=Nsnap;
	Cat->Ngroups=SubCat.Ngroups;
	Cat->Nsubs=SubCat.Nsubs;
	Cat->SubHalo=mymalloc(sizeof(SubData)*Cat->Nsubs);
	if(Nsnap<MaxSnap-1)
	load_pro2dest(Nsnap,&pro2dest,&Npro,SUBCAT_DIR);
	else  //forge pro2dest for the last snapshot
	{
		pro2dest=mymalloc(sizeof(HBTInt)*(SubCat.Nsubs+1));
		pro2dest++;
		for(i=0;i<SubCat.Nsubs;i++)
			pro2dest[i]=-2;//unknown descendent
		pro2dest[-1]=-1;
		Npro=SubCat.Nsubs;
	}
	if(Npro<Cat->Nsubs)	{printf("pro2dest erro!\n");exit(1);}
	for(i=0;i<Cat->Nsubs;i++)
	{
	  Cat->SubHalo[i].Mdm=SubCat.SubLen[i];
	  Cat->SubHalo[i].HostID=SubCat.HaloChains[i].HostID;
	  Cat->SubHalo[i].ProID=SubCat.HaloChains[i].ProSubID;
	  Cat->SubHalo[i].DesID=pro2dest[i];
	  Cat->SubHalo[i].SubRank=SubCat.SubRank[i];
	}
	free_pro2dest(pro2dest);
	free_sub_table(&SubCat);
	
	return sizeof(SubData)*Cat->Nsubs+3*sizeof(HBTInt);
}
void erase_brfcat(BRFCAT *Cat)
{
	myfree(Cat->SubHalo);
}

void save_sub2hist(HBTInt Nsnap,HBTInt *Sub2Hist,HBTInt Nsubs,char *subcatdir)
{
FILE *fp;
char buf[1024];

  sprintf(buf, "%s/history", subcatdir);
  mkdir(buf,0755);
  sprintf(buf, "%s/history/sub2hist_%03d", subcatdir, (int)Nsnap);
  myfopen(fp,buf,"w");
  fwrite(&Nsubs,sizeof(HBTInt),1,fp);
  fwrite(Sub2Hist,sizeof(HBTInt),Nsubs,fp);
  fclose(fp);	
}
void load_sub2hist(HBTInt Nsnap,HBTInt **P2sub2hist,HBTInt *P2Nsubs,char *subcatdir)
{
FILE *fp;
char buf[1024];
HBTInt Nsubs,*Sub2Hist;

  sprintf(buf, "%s/history/sub2hist_%03d", subcatdir, (int)Nsnap);
  myfopen(fp,buf,"r");
  fread(&Nsubs,sizeof(HBTInt),1,fp);
  Sub2Hist=mymalloc(sizeof(HBTInt)*Nsubs);
  fread(Sub2Hist,sizeof(HBTInt),Nsubs,fp);
  fclose(fp);
  *P2Nsubs=Nsubs;
  *P2sub2hist=Sub2Hist;
}
void save_evocat_pre(EVOLUTIONCAT_Pre *EvoCat,char *subcatdir)
{
FILE *fp;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/history", subcatdir);
  mkdir(buf,0755);
  sprintf(buf, "%s/history/EvoCat_pre", subcatdir);
  myfopen(fp,buf,"w");
  fwrite(&EvoCat->NNode,sizeof(HBTInt),1,fp);
  fwrite(&EvoCat->NHist,sizeof(HBTInt),1,fp);
  fwrite(EvoCat->HistLen,sizeof(HBTInt),EvoCat->NHist,fp);
  fwrite(EvoCat->HistOffset,sizeof(HBTInt),EvoCat->NHist,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].SnapBirth,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].SnapDeath,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].SnapEnter,sizeof(HBTInt),1,fp); 
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].ProHistID,sizeof(HBTInt),1,fp); 	 
 fwrite(EvoCat->History[0].Member,sizeof(SubNodePre),EvoCat->NNode,fp);	 
  fclose(fp);		
}
void load_evocat_pre(EVOLUTIONCAT_Pre *EvoCat,char *subcatdir)
{
FILE *fp;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/history/EvoCat_pre", subcatdir);
  myfopen(fp,buf,"r");
  fread(&EvoCat->NNode,sizeof(HBTInt),1,fp);
  fread(&EvoCat->NHist,sizeof(HBTInt),1,fp);
  EvoCat->HistLen=mymalloc(sizeof(HBTInt)*EvoCat->NHist);
  EvoCat->HistOffset=mymalloc(sizeof(HBTInt)*EvoCat->NHist);
  EvoCat->History=mymalloc(sizeof(HISTORY_Pre)*EvoCat->NHist);
  EvoCat->History[0].Member=mymalloc(sizeof(SubNode)*EvoCat->NNode);
  fread(EvoCat->HistLen,sizeof(HBTInt),EvoCat->NHist,fp);
  fread(EvoCat->HistOffset,sizeof(HBTInt),EvoCat->NHist,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapBirth,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapDeath,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapEnter,sizeof(HBTInt),1,fp); 
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].ProHistID,sizeof(HBTInt),1,fp); 	 
 	
 fread(EvoCat->History[0].Member,sizeof(SubNode),EvoCat->NNode,fp);	 
 for(i=0;i<EvoCat->NHist;i++)
 	EvoCat->History[i].Member=EvoCat->History[0].Member+EvoCat->HistOffset[i];
  fclose(fp);
}
void free_evocat_pre(EVOLUTIONCAT_Pre *ECat)
{
	myfree(ECat->HistLen);
	myfree(ECat->HistOffset);
	myfree(ECat->History[0].Member);
	myfree(ECat->History);
}
void load_history_pre(EVOLUTIONCAT *EvoCat,char *subcatdir)
{//same as load_evocat_pre, except for using EVOLUTIONCAT rather than EVOLUTIONCAT_Pre
FILE *fp;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/history/EvoCat_pre", subcatdir);
  myfopen(fp,buf,"r");
  fread(&EvoCat->NNode,sizeof(HBTInt),1,fp);
  fread(&EvoCat->NHist,sizeof(HBTInt),1,fp);
  EvoCat->HistLen=mymalloc(sizeof(HBTInt)*EvoCat->NHist);
  EvoCat->HistOffset=mymalloc(sizeof(HBTInt)*EvoCat->NHist);
  EvoCat->History=mymalloc(sizeof(HISTORY)*EvoCat->NHist);
  EvoCat->History[0].Member=mymalloc(sizeof(SubNode)*EvoCat->NNode);
  fread(EvoCat->HistLen,sizeof(HBTInt),EvoCat->NHist,fp);
  fread(EvoCat->HistOffset,sizeof(HBTInt),EvoCat->NHist,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapBirth,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapDeath,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapEnter,sizeof(HBTInt),1,fp); 
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].ProHistID,sizeof(HBTInt),1,fp); 	 
 for(i=0;i<EvoCat->NNode;i++)
 {	
 	fread(&EvoCat->History[0].Member[i].Mdm,sizeof(HBTInt),1,fp);	 
	fread(&EvoCat->History[0].Member[i].SubID,sizeof(HBTInt),1,fp);	 
	fread(&EvoCat->History[0].Member[i].HostID,sizeof(HBTInt),1,fp);	 
	fread(&EvoCat->History[0].Member[i].SubRank,sizeof(HBTInt),1,fp);	 
 }
 for(i=0;i<EvoCat->NHist;i++)
 	EvoCat->History[i].Member=EvoCat->History[0].Member+EvoCat->HistOffset[i];
  fclose(fp);
}

void save_evocat_rev(EVOLUTIONCAT *EvoCat,char *subcatdir)
{
FILE *fp;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/history", subcatdir);
  mkdir(buf,0755);
  sprintf(buf, "%s/history/EvoCat_rev", subcatdir);
  myfopen(fp,buf,"w");
  fwrite(&EvoCat->PartMass,sizeof(HBTReal),1,fp);
  fwrite(&EvoCat->NNode,sizeof(HBTInt),1,fp);
  fwrite(&EvoCat->NHist,sizeof(HBTInt),1,fp);
  fwrite(EvoCat->HistLen,sizeof(HBTInt),EvoCat->NHist,fp);
  fwrite(EvoCat->HistOffset,sizeof(HBTInt),EvoCat->NHist,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].SnapBirth,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].SnapDeath,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].SnapEnter,sizeof(HBTInt),1,fp); 
  for(i=0;i<EvoCat->NHist;i++)
 	 fwrite(&EvoCat->History[i].ProHistID,sizeof(HBTInt),1,fp); 	 
 fwrite(EvoCat->History[0].Member,sizeof(SubNode),EvoCat->NNode,fp);	 
  fclose(fp);		
}
void load_evocat_rev(EVOLUTIONCAT *EvoCat,char *subcatdir)
{
FILE *fp;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/history/EvoCat_rev", subcatdir);
  myfopen(fp,buf,"r");
  fread(&EvoCat->PartMass,sizeof(HBTReal),1,fp);
  fread(&EvoCat->NNode,sizeof(HBTInt),1,fp);
  fread(&EvoCat->NHist,sizeof(HBTInt),1,fp);
  EvoCat->HistLen=mymalloc(sizeof(HBTInt)*EvoCat->NHist);
  EvoCat->HistOffset=mymalloc(sizeof(HBTInt)*EvoCat->NHist);
  EvoCat->History=mymalloc(sizeof(HISTORY)*EvoCat->NHist);
  EvoCat->History[0].Member=mymalloc(sizeof(SubNode)*EvoCat->NNode);
  fread(EvoCat->HistLen,sizeof(HBTInt),EvoCat->NHist,fp);
  fread(EvoCat->HistOffset,sizeof(HBTInt),EvoCat->NHist,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapBirth,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapDeath,sizeof(HBTInt),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapEnter,sizeof(HBTInt),1,fp); 
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].ProHistID,sizeof(HBTInt),1,fp); 	 
 	
 fread(EvoCat->History[0].Member,sizeof(SubNode),EvoCat->NNode,fp);	 
 for(i=0;i<EvoCat->NHist;i++)
 	EvoCat->History[i].Member=EvoCat->History[0].Member+EvoCat->HistOffset[i];
  fclose(fp);
}

SubNode *GetMember(EVOLUTIONCAT *EvoCat,HBTInt HistID,HBTInt Nsnap)
{
	HBTInt SnapBirth,SnapDeath,ProHistID;
	SnapBirth=EvoCat->History[HistID].SnapBirth;
	SnapDeath=EvoCat->History[HistID].SnapDeath;

	if(Nsnap<SnapBirth)
		if((ProHistID=EvoCat->History[HistID].ProHistID)>0)//splitter,have members before birth
			return GetMember(EvoCat,ProHistID,Nsnap);
		else
			return NULL;

	if(Nsnap>=SnapDeath)
		return NULL;

	return &(EvoCat->History[HistID].Member[Nsnap-SnapBirth]);
}

void create_historyshards(struct HistoryShards *HistoryShard,HBTInt NumHist)
{
HistoryShard->NumHist=NumHist;
HistoryShard->NBirth=mymalloc(sizeof(HBTInt)*NumHist);
HistoryShard->HistoryOffset=mymalloc(sizeof(HBTInt)*NumHist);
HistoryShard->Par=mymalloc(sizeof(struct ShardParam *)*NumHist);
}
void save_historyshards(struct HistoryShards *HistoryShard)
{
HBTInt NumHist,i,Offset;
char buf[1024];
FILE *fp;
  
  NumHist=HistoryShard->NumHist;
  Offset=0;	
  for(i=0;i<NumHist;i++)
  {
  HistoryShard->HistoryOffset[i]=Offset;
  Offset+=HistoryShard->NBirth[i];
  }
  HistoryShard->NumShards=Offset;
  	
  sprintf(buf, "%s/history/HistoryShards", SUBCAT_DIR);
  myfopen(fp,buf,"w");
  fwrite(&NumHist,sizeof(HBTInt),1,fp);
  fwrite(&Offset,sizeof(HBTInt),1,fp);
  fwrite(HistoryShard->NBirth,sizeof(HBTInt),NumHist,fp);
  fwrite(HistoryShard->HistoryOffset,sizeof(HBTInt),NumHist,fp);
  for(i=0;i<NumHist;i++)
  fwrite(HistoryShard->Par[i],sizeof(struct ShardParam),HistoryShard->NBirth[i],fp);
  fclose(fp);	
}	
void load_historyshards(struct HistoryShards *HistoryShard)
{
HBTInt NumHist,i;
char buf[1024];
FILE *fp;
  	
  sprintf(buf, "%s/history/HistoryShards", SUBCAT_DIR);
  myfopen(fp,buf,"r");
  fread(&NumHist,sizeof(HBTInt),1,fp);
  
  HistoryShard->NumHist=NumHist;
  HistoryShard->NBirth=mymalloc(sizeof(HBTInt)*NumHist);
  HistoryShard->HistoryOffset=mymalloc(sizeof(HBTInt)*NumHist);
  HistoryShard->Par=mymalloc(sizeof(struct ShardParam *)*NumHist);
  
  fread(&HistoryShard->NumShards,sizeof(HBTInt),1,fp);
  fread(HistoryShard->NBirth,sizeof(HBTInt),NumHist,fp);
  fread(HistoryShard->HistoryOffset,sizeof(HBTInt),NumHist,fp);
  for(i=0;i<NumHist;i++)
  {
	HistoryShard->Par[i]=mymalloc(sizeof(struct ShardParam)*HistoryShard->NBirth[i]);  
	fread(HistoryShard->Par[i],sizeof(struct ShardParam),HistoryShard->NBirth[i],fp);
  }
  fclose(fp);	
}			

void fresh_MBDID2Index(MBDCATALOGUE *MbdCat)
{
	HBTInt i,pid;			
	
	if(0==PIDHash.np)
	{
		fprintf(logfile,"call fill_PIDHash() before fresh_ID2Index()!\n");
		exit(1);
	}

	for(i=0;i<MbdCat->NSubs+MbdCat->NOrphans;i++)
	{
		pid=MbdCat->Nodes[i].MBD_PID;
		MbdCat->Nodes[i].MBD_PID=lookup_ID2Ind(pid);
	}
}

void save_mbd_catalogue(HBTInt Nsnap, MBDCATALOGUE *MbdCat)
{
FILE *fd;
char buf[1024];
HBTInt i,j;

  fprintf(logfile,"saving most-bound particles..\n");
  fprintf(logfile,"Nsub="HBTIFMT", NOrphan="HBTIFMT", NQuasi="HBTIFMT", NOrphanQuasi="HBTIFMT"\n",
				MbdCat->NSubs,MbdCat->NOrphans,MbdCat->NQuasi,MbdCat->NOrphanQuasi);
  fflush(logfile);
  sprintf(buf, "%s/mostbound/mbdcat_%03d", SUBCAT_DIR, (int)Nsnap);
  myfopen(fd,buf,"w");

  fwrite(&MbdCat->NSubs,sizeof(HBTInt),1,fd);
  fwrite(&MbdCat->NOrphans,sizeof(HBTInt),1,fd);
  fwrite(&MbdCat->Ngroups,sizeof(HBTInt),1,fd);
  fwrite(&MbdCat->NQuasi,sizeof(HBTInt),1,fd);
  fwrite(&MbdCat->NOrphanQuasi,sizeof(HBTInt),1,fd);
  fwrite(MbdCat->GrpLen_Sub,sizeof(HBTInt),MbdCat->Ngroups,fd);
  fwrite(MbdCat->GrpOffset_Sub,sizeof(HBTInt),MbdCat->Ngroups,fd);
  fwrite(MbdCat->GrpLen_Orphan,sizeof(HBTInt),MbdCat->Ngroups,fd);
  fwrite(MbdCat->GrpOffset_Orphan,sizeof(HBTInt),MbdCat->Ngroups,fd);
  for(i=0;i<MbdCat->NSubs+MbdCat->NOrphans;i++)
  {
	fwrite(&(MbdCat->Nodes[i].MBD_PID),sizeof(HBTInt),1,fd);
	fwrite(&(MbdCat->Nodes[i].HistID),sizeof(HBTInt),1,fd);
	fwrite(&(MbdCat->Nodes[i].HostID),sizeof(HBTInt),1,fd);
	fwrite(MbdCat->Nodes[i].Pos,sizeof(HBTReal),3,fd);
	fwrite(MbdCat->Nodes[i].Vel,sizeof(HBTReal),3,fd);
  }

fclose(fd);	
}

void load_mbd_catalogue(HBTInt Nsnap, MBDCATALOGUE *MbdCat)
{
FILE *fd;
char buf[1024];
HBTInt i,j;

  sprintf(buf, "%s/mostbound/mbdcat_%03d", SUBCAT_DIR, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&MbdCat->NSubs,sizeof(HBTInt),1,fd);
  fread(&MbdCat->NOrphans,sizeof(HBTInt),1,fd);
  fread(&MbdCat->Ngroups,sizeof(HBTInt),1,fd);
  fread(&MbdCat->NQuasi,sizeof(HBTInt),1,fd);
  fread(&MbdCat->NOrphanQuasi,sizeof(HBTInt),1,fd);
  
  MbdCat->GrpLen_Sub=mymalloc(sizeof(HBTInt)*MbdCat->Ngroups);
  MbdCat->GrpOffset_Sub=mymalloc(sizeof(HBTInt)*MbdCat->Ngroups);
  MbdCat->GrpLen_Orphan=mymalloc(sizeof(HBTInt)*MbdCat->Ngroups);
  MbdCat->GrpOffset_Orphan=mymalloc(sizeof(HBTInt)*MbdCat->Ngroups);
MbdCat->Nodes=mymalloc(sizeof(MBDNode)*(MbdCat->NSubs+MbdCat->NOrphans));

  fread(MbdCat->GrpLen_Sub,sizeof(HBTInt),MbdCat->Ngroups,fd);
  fread(MbdCat->GrpOffset_Sub,sizeof(HBTInt),MbdCat->Ngroups,fd);
  fread(MbdCat->GrpLen_Orphan,sizeof(HBTInt),MbdCat->Ngroups,fd);
  fread(MbdCat->GrpOffset_Orphan,sizeof(HBTInt),MbdCat->Ngroups,fd);
  
  for(i=0;i<MbdCat->NSubs+MbdCat->NOrphans;i++)
  {
	fread(&MbdCat->Nodes[i].MBD_PID,sizeof(HBTInt),1,fd);
	fread(&MbdCat->Nodes[i].HistID,sizeof(HBTInt),1,fd);
	fread(&MbdCat->Nodes[i].HostID,sizeof(HBTInt),1,fd);
	fread(MbdCat->Nodes[i].Pos,sizeof(HBTReal),3,fd);
	fread(MbdCat->Nodes[i].Vel,sizeof(HBTReal),3,fd);
  }

fclose(fd);	
}

void free_mbd_catalogue(MBDCATALOGUE *MbdCat)
{
	myfree(MbdCat->GrpLen_Sub);
	myfree(MbdCat->GrpOffset_Sub);
	myfree(MbdCat->GrpLen_Orphan);
	myfree(MbdCat->GrpOffset_Orphan);
	myfree(MbdCat->Nodes);
}
