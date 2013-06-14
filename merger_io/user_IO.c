#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define myfread(a,b,c,d) fread_swap(a,b,c,d,ByteOrder)

#ifdef SKIP
	#undef SKIP
	#undef SKIP2
	#undef CHECK
#endif
#define SKIP myfread(&dummy,sizeof(dummy),1,fp)
#define SKIP2 myfread(&dummy2,sizeof(dummy2),1,fp)
#define CHECK if(dummy!=dummy2){fprintf(logfile,"error!record brackets not match for Snap%d!\t%d,%d\n",(int)Nsnap,dummy,dummy2);fflush(logfile);exit(1);} 

int check_snapshot_byteorder(char *filename)
{
	/* to check whether byteswap is needed, return 1 if yes, 0 if no, exit if error*/
	int dummy,n,ns;
	FILE *fp;

	n=sizeof(header);
	ns=n;
	swap_Nbyte(&ns,1,sizeof(ns));
	myfopen(fp,filename,"r");
	fread(&dummy,sizeof(dummy),1,fp);
	fclose(fp);
	
	if(dummy==n) return 0;
	if(dummy==ns) return 1;
	
	fprintf(logfile,"endianness check failed for: %s \n file format not expected:%d;%d,%d\n", filename,dummy,n,ns);
	fflush(logfile);
	exit(1);
}

int check_grpcat_byteorder(char *filename)
{
	/* to check whether byteswap is needed, return 1 if yes, 0 if no, exit if error*/
	int Nfiles,n,ns;
	long offset;
	FILE *fp;

	n=NFILES_GRP;
	ns=n;
	swap_Nbyte(&ns,1,sizeof(ns));
	
	myfopen(fp,filename,"r");
	#ifdef GRP_V3FORMAT
	offset=5*4;  //3*int+1*longlong
	#else
	offset=3*4;  //3*int
	#endif
	fseek(fp,offset,SEEK_SET);
	fread(&Nfiles,sizeof(int),1,fp);
	fclose(fp);
	
	if(Nfiles==n) return 0;
	if(Nfiles==ns) return 1;

	fprintf(logfile,"endianness check failed for: %s , file format not expected:%d;%d,%d\n", filename,(int)Nfiles,n,ns);
	fflush(logfile);
	exit(1);
}

static int read_gadget_header(FILE *fp,IO_HEADER * h)
{//read the header part, assign header extensions, and do several consistency check
//return ByteOrder
	int dummy,dummy2,n,ns,ByteOrder;
	
	n=sizeof(IO_HEADER);
	ns=n;
	swap_Nbyte(&ns,1,sizeof(ns));
			
	fread(&dummy,sizeof(dummy),1,fp);
	
	if(dummy==n)
	 ByteOrder=0;
	else if(dummy==ns)
	 ByteOrder=1;
	else
	{
		fprintf(logfile,"endianness check failed for header\n file format not expected:%d;%d,%d\n",dummy,n,ns);
		fflush(logfile);
		exit(1);
	}
	
	dummy=n;

	myfread(h->npart,sizeof(int),6,fp);
	myfread(h->mass,sizeof(double),6,fp);
	myfread(&h->time,sizeof(double),1,fp);
	myfread(&h->redshift,sizeof(double),1,fp);
	myfread(&h->flag_sfr,sizeof(int),1,fp);
	myfread(&h->flag_feedback,sizeof(int),1,fp);
	myfread(h->npartTotal,sizeof(int),6,fp);
	myfread(&h->flag_cooling,sizeof(int),1,fp);
	myfread(&h->num_files,sizeof(int),1,fp);
	myfread(&h->BoxSize,sizeof(double),1,fp);
	myfread(&h->Omega0,sizeof(double),1,fp);
	myfread(&h->OmegaLambda,sizeof(double),1,fp);
	myfread(&h->HubbleParam,sizeof(double),1,fp);
	fseek(fp,n+sizeof(int),SEEK_SET);
	myfread(&dummy2,sizeof(dummy2),1,fp);
    if(dummy!=dummy2)
	{
		fprintf(logfile,"error!record brackets not match for header!\t%d,%d\n",dummy,dummy2);
		exit(1);
	} 
	
	/*extend and examine the header*/
	if(NFILES!=h->num_files)
	  {
		  fprintf(logfile,"error: number of snapfiles specified not the same as stored: %d,%d\n",
		  NFILES,h->num_files);
		  fflush(logfile);
//		  exit(1);
	  }
	  
	h->Hz=HUBBLE0 * sqrt(h->Omega0 / (h->time * h->time * h->time) 
			+ (1 - h->Omega0 - h->OmegaLambda) / (h->time * h->time)
			+ h->OmegaLambda);//Hubble param for the current catalogue;
	
	if(h->npartTotal[0])
	  {
	  if(0==h->mass[0])
	  h->mass[0]=MP_GAS;
	  else if(fabs(MP_GAS-h->mass[0])>0.01*MP_GAS)
	  {
		  fprintf(logfile,"error: GAS particle mass mismatch: %g,%g\n",MP_GAS,h->mass[0]);
		  fflush(logfile);
		  exit(1);
	  }
	  }
     if(0==h->mass[1])
	  h->mass[1]=MP_DM;
	  else if(fabs(MP_DM-h->mass[1])>0.01*MP_DM)
	  {
		  fprintf(logfile,"error: DM particle mass mismatch: %g,%g\n",MP_DM,h->mass[1]);
		  fflush(logfile);
		  exit(1);
	  }
	  
	#ifdef CONVERT_LENGTH_MPC_KPC
	  h->BoxSize*=1000.;
	#endif 
	
	int j;
	for(j=0, NumPart=0; j<6; j++)
	    NumPart+= h->npartTotal[j];	
	if(NP_SIM!=NumPart)
	{
		fprintf(logfile,"error:Total Number of Particles differ: %lld,%lld\n",(long long)NP_SIM,(long long)NumPart);
		fflush(logfile);
		exit(1);
	}
	if(NP_GAS!=h->npartTotal[0])
	{
		fprintf(logfile,"error:Number of Gas Particles differ: %lld,%u\n",(long long)NP_GAS,h->npartTotal[0]);
		fflush(logfile);
		exit(1);
	}
	if(NP_DM!=h->npartTotal[1])
	{
		fprintf(logfile,"error:Number of DM Particles differ: %lld,%u\n",(long long)NP_DM,h->npartTotal[1]);
		fflush(logfile);
		exit(1);
	}
	
	return ByteOrder;
	
}
void load_particle_header(HBTInt Nsnap, char *SnapPath)
{
	FILE *fp;
	char buf[1024];
	HBTInt NumPart;

	#ifdef SNAPLIST
	Nsnap=snaplist[Nsnap];
	#endif
	sprintf(buf,"%s/snapdir_%03d/%s_%03d.0",SnapPath,(int)Nsnap,SNAPFILE_BASE,(int)Nsnap);
    if(1==NFILES)
	 if(!try_readfile(buf))	sprintf(buf,"%s/%s_%03d",SnapPath,SNAPFILE_BASE,(int)Nsnap); //try the other convention
	
	myfopen(fp,buf,"r");
	read_gadget_header(fp,&header);
	fclose(fp);
	header.Nsnap=Nsnap;
}

struct io_ID2Ind
  {
	IDatInt PID;
	HBTInt PInd;  //address in Pdat.PID
  };
static int comp_PIDArr(const void *a, const void *b)//used to sort PID in ascending order; 
{
  if(((struct io_ID2Ind *) a)->PID > ((struct io_ID2Ind *) b)->PID)
    return +1;

  if(((struct io_ID2Ind *) a)->PID < ((struct io_ID2Ind *) b)->PID)
    return -1;

  return 0;
}
static int comp_IDatInt(const void *a, const void *b)//used to sort PID in ascending order; 
{
  if((*(IDatInt *)a) > (*(IDatInt *)b))
    return +1;

  if((*(IDatInt *)a) < (*(IDatInt *)b))
    return -1;

  return 0;
}
/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
void load_particle_data(HBTInt Nsnap, char *SnapPath)
{
	int dummy,dummy2,ByteOrder;
	long int pre_len,tail_len;
	FILE *fp;
	char buf[1024];
	HBTInt i,j;
	HBTInt Nload=0;
	struct 
	{
	IDatInt *PID;	
	IDatReal (*Pos)[3];
	IDatReal (*Vel)[3];
	} IDat; //input particle data, temporary.
	
	IDat.PID=mymalloc(sizeof(IDatInt)*NP_DM);
	IDat.Pos=mymalloc(sizeof(IDatReal)*3*NP_DM);
	IDat.Vel=mymalloc(sizeof(IDatReal)*3*NP_DM);
	
	Pdat.Nsnap=Nsnap;
	#ifdef SNAPLIST
	Nsnap=snaplist[Nsnap];
	#endif

  for(i=0; i<NFILES; i++)
    {
	sprintf(buf,"%s/snapdir_%03d/%s_%03d.%d",SnapPath,(int)Nsnap,SNAPFILE_BASE,(int)Nsnap,(int)i);
    if(1==NFILES)
	 if(!try_readfile(buf))	sprintf(buf,"%s/%s_%03d",SnapPath,SNAPFILE_BASE,(int)Nsnap); //try the other convention

      myfopen(fp,buf,"r");

      ByteOrder=read_gadget_header(fp,&header);
	  header.Nsnap=Nsnap;
	  
	  pre_len=header.npart[0]*sizeof(IDatReal)*3;
	  for(j=2,tail_len=0;j<6;j++)tail_len+=header.npart[j];
	  tail_len*=sizeof(IDatReal)*3;
      SKIP;
	  fseek(fp,pre_len,SEEK_CUR);//load only DM data;
	  myfread(IDat.Pos+Nload,sizeof(IDatReal)*3,header.npart[1],fp);
	  fseek(fp,tail_len,SEEK_CUR);
      SKIP2;
	  CHECK;
	
      SKIP;
      fseek(fp,pre_len,SEEK_CUR);//load only DM data;
	  myfread(IDat.Vel+Nload,sizeof(IDatReal)*3,header.npart[1],fp);
	  fseek(fp,tail_len,SEEK_CUR);
      SKIP2;
	  CHECK;
      
	  pre_len=header.npart[0]*sizeof(IDatInt);
	  for(j=2,tail_len=0;j<6;j++)tail_len+=header.npart[j];
	  tail_len*=sizeof(IDatInt);
	  SKIP;
      fseek(fp,pre_len,SEEK_CUR);//load only DM data;
	  myfread(IDat.PID+Nload,sizeof(IDatInt),header.npart[1],fp);
	  fseek(fp,tail_len,SEEK_CUR);
      SKIP2;
	  CHECK;
	  
	  if(feof(fp))
	  {
		fprintf(logfile,"error:End-of-File in %s\n",buf);
		fflush(logfile);exit(1);  
	  }
	
	  Nload+=header.npart[1];
     
      fclose(fp);
  }
  if(NP_DM!=Nload)
  {
	  fprintf(logfile,"error: Number of loaded DM particles mismatch: %lld,%lld\n",(long long)NP_DM,(long long)Nload);
	  fflush(logfile);
	  exit(1);
  }
  
//now transfer to HBT's internal data  
  #ifdef HBTPID_RANKSTYLE
  struct io_ID2Ind *table;
  table=mymalloc(sizeof(struct io_ID2Ind)*NP_DM);
  for(i=0;i<NP_DM;i++)
  {
	  table[i].PID=IDat.PID[i];
	  table[i].PInd=i;
  }
  qsort(table,NP_DM,sizeof(struct io_ID2Ind),comp_PIDArr);
  Pdat.PID=mymalloc(sizeof(HBTInt)*NP_DM);
  for(i=0;i<NP_DM;i++)
	Pdat.PID[table[i].PInd]=i;  //now PID has been turned into particle ranks
	
  sprintf(buf,"%s/DM_PIDs_Sorted.dat",SUBCAT_DIR);
  if(!try_readfile(buf))//create the file if it does not exist
  {
	  IDatInt np;
	  myfopen(fp,buf,"w");
	  np=NP_DM;
	  fwrite(&np,sizeof(IDatInt),1,fp);
	  for(i=0;i<NP_DM;i++)
	  fwrite(&(table[i].PID),sizeof(IDatInt),1,fp);
	  fwrite(&np,sizeof(IDatInt),1,fp);
	  fclose(fp);
  }
	
  myfree(table);  
  myfree(IDat.PID);
  #else
#ifdef SAME_INTTYPE
  Pdat.PID=IDat.PID;
#else
  Pdat.PID=mymalloc(sizeof(HBTInt)*NP_DM);
  for(i=0;i<NP_DM;i++)
	Pdat.PID[i]=IDat.PID[i];
  myfree(IDat.PID);
#endif
  #endif
#ifdef SAME_REALTYPE
	  Pdat.Pos=IDat.Pos;
	  Pdat.Vel=IDat.Vel;
#else
	  Pdat.Pos=mymalloc(sizeof(HBTxyz)*NP_DM);
	  Pdat.Vel=mymalloc(sizeof(HBTxyz)*NP_DM);
	  for(i=0;i<NP_DM;i++)
	  {
		for(j=0;j<3;j++)
		{
			Pdat.Pos[i][j]=IDat.Pos[i][j];
			Pdat.Vel[i][j]=IDat.Vel[i][j];
		}
      }
	  myfree(IDat.Pos);
	  myfree(IDat.Vel);
#endif
  	
  #ifdef CONVERT_LENGTH_MPC_KPC
  for(i=0;i<NP_DM;i++)
  for(j=0;j<3;j++)
  Pdat.Pos[i][j]*=1000.;
  #endif	
}

void free_particle_data()
{
#ifndef PID_ORDERED
	myfree(Pdat.PID);
#endif
	myfree(Pdat.Pos);
	myfree(Pdat.Vel);
	Pdat.Nsnap=-100;
}

IDatInt *load_PIDs_Sorted()
{
#ifdef HBTPID_RANKSTYLE
  char buf[1024];
  FILE *fd;
  IDatInt *PIDs,np,np2;  
  sprintf(buf,"%s/DM_PIDs_Sorted.dat",SUBCAT_DIR);
  if(!try_readfile(buf))//create the file if it does not exist
  {
	 load_particle_data(0,SNAPSHOT_DIR);
	 free_particle_data();
  }
  myfopen(fd,buf,"r");
  fread(&np,sizeof(IDatInt),1,fd);
  if(np!=NP_DM){fprintf(logfile,"error, file corruption: %s\n",buf);exit(1);}
  
  PIDs=mymalloc(sizeof(IDatInt)*NP_DM);
  fread(PIDs,sizeof(IDatInt),NP_DM,fd);
  
  fread(&np2,sizeof(IDatInt),1,fd);
  if(np2!=NP_DM){fprintf(logfile,"error, file corruption: %s\n",buf);exit(1);}
  fclose(fd);

  return PIDs;
#else
  return NULL;
#endif
}

#ifdef LOAD_SNAP_AS_GRP  //intended for major merger tests.
void load_group_catalogue(HBTInt Nsnap,CATALOGUE *Cat, char *GrpPath)
{
	fprintf(logfile,"Nsnap=%d\n",Nsnap);
	if(Pdat.Nsnap!=Nsnap)
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	
	if(0==Nsnap)
	{
	  Cat->Ngroups=2;
	  Cat->Nids=NP_SIM;
	  Cat->Len=mymalloc(sizeof(HBTInt)*2);
	  Cat->Offset=mymalloc(sizeof(HBTInt)*2);
	  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*NP_DM);
	  Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);
	  HBTInt i,len0, len1, *PID2;
	  PID2=mymalloc(sizeof(HBTInt)*NP_DM);
          for(i=0,len0=0,len1=0;i<NP_DM;i++)
          {
                if(Pdat.Pos[i][0]<55.)
                {
                Cat->PIDorIndex[len0]=Pdat.PID[i];
                len0++;
                }
                else
                {
                PID2[len1]=Pdat.PID[i];
                len1++;
                }
          }
	Cat->Len[0]=len0;
	Cat->Len[1]=len1;
	if(NP_DM!=len0+len1)
	{
	printf("error in data division: %d + %d != %d\n", len0, len1, NP_DM);
	exit(1);
	}
	memcpy(Cat->PIDorIndex+len0,PID2,sizeof(HBTInt)*len1);
	myfree(PID2);
	Cat->Offset[0]=0;
	Cat->Offset[1]=len0;
	}
	else
	{
	Cat->Ngroups=1;
	Cat->Nids=NP_SIM;
	Cat->Len=mymalloc(sizeof(HBTInt));
	Cat->Offset=mymalloc(sizeof(HBTInt));
	Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*NP_DM);
	Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);
	Cat->Len[0]=NP_SIM;
	Cat->Offset[0]=0;
	memcpy(Cat->PIDorIndex,Pdat.PID,sizeof(HBTInt)*NP_DM);
	}
}
#else
void load_group_catalogue(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{//wrapper for v2 or v3
#ifdef GRP_V3FORMAT
load_group_catalogue_v3(Nsnap,Cat,GrpPath);
#else
	#ifdef GRP_MILLIFORMAT
load_group_catalogue_millimill(Nsnap,Cat,GrpPath);
	#else
load_group_catalogue_v2(Nsnap,Cat,GrpPath);	
	#endif
#endif
}
#endif
	
void load_group_catalogue_v3(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{//PGADGET-3's subfind format
  FILE *fd;
  char buf[1024];
  int Ngroups,TotNgroups,Nids,NFiles;
  long long i,Nload,TotNids;
  int ByteOrder,grpfile_type;
  
struct 
{
int *Len;
int *Offset;
IDatInt *PID; 
}ICat;
  
  	#ifdef SNAPLIST
	Nsnap=snaplist[Nsnap];
	#endif

  Nload=0;
  for(i=0;i<NFILES_GRP;i++)
  {
  grpfile_type=1;	
  sprintf(buf, "%s/groups_%03d/subhalo_tab_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
  if(!try_readfile(buf))
  {
	  sprintf(buf, "%s/groups_%03d/group_tab_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	  grpfile_type=2;
  }
  if(1==NFILES_GRP)
  {
	if(!try_readfile(buf))
	{
		sprintf(buf, "%s/subhalo_tab_%03d",GrpPath,(int)Nsnap);
		grpfile_type=3;
	}
	if(!try_readfile(buf))
	{
	sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
	grpfile_type=4;
	}
  }
	
  ByteOrder=check_grpcat_byteorder(buf);
  	
  myfopen(fd,buf,"r");

  myfread(&Ngroups, sizeof(int), 1, fd);
  myfread(&TotNgroups, sizeof(int), 1, fd);
  Cat->Ngroups=TotNgroups;
  myfread(&Nids, sizeof(int), 1, fd);
  myfread(&TotNids,sizeof(long long),1,fd);
  #ifndef HBT_INT8
  if(TotNids>INT_MAX)
  {
	  printf("error: TotNids larger than HBTInt can hold! %lld\n",TotNids);
	  exit(1);
  }
  #endif
  Cat->Nids=TotNids;
  myfread(&NFiles, sizeof(int), 1, fd);
 
  if(1==grpfile_type||3==grpfile_type)
  {
  int Nsub,TotNsub;
  myfread(&Nsub,sizeof(int),1,fd);
  myfread(&TotNsub,sizeof(int),1,fd);
  }
  
  if(NFILES_GRP!=NFiles)
	  {
		  fprintf(logfile,"error: number of grpfiles specified not the same as stored: %d,%d\n for file %s\n",
		  NFILES_GRP,(int)NFiles,buf);
		  fflush(logfile);
		  exit(1);
	  }
 
  if(0==i)
  fprintf(logfile,"Snap="HBTIFMT"  Ngroups=%d  Nids=%d  TotNgroups=%d  NFiles=%d\n", Nsnap, Ngroups, Nids, (int)(Cat->Ngroups), NFiles);
  
  if(0==i)
  {
  ICat.Len= mymalloc(sizeof(int)*Cat->Ngroups);
  ICat.Offset=mymalloc(sizeof(int)*Cat->Ngroups);
  }
  
  myfread(ICat.Len+Nload, sizeof(int), Ngroups, fd);
  myfread(ICat.Offset+Nload, sizeof(int), Ngroups, fd);
  if(feof(fd))
  {
	fprintf(logfile,"error:End-of-File in %s\n",buf);
	fflush(logfile);exit(1);  
  }
  Nload+=Ngroups;
  fclose(fd);
  }
  if(Nload!=Cat->Ngroups)
  {
	  fprintf(logfile,"error:Num groups loaded not match: %lld,"HBTIFMT"\n",Nload,Cat->Ngroups);
	  fflush(logfile);exit(1);
  }
  
  ICat.PID=mymalloc(sizeof(IDatInt)*Cat->Nids);
  Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);//consider move this out.............
  Nload=0;
  for(i=0;i<NFILES_GRP;i++)
  {
  switch(grpfile_type)
  {
  case 1:
	sprintf(buf, "%s/groups_%03d/subhalo_ids_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	break;
  case 2:
	sprintf(buf, "%s/groups_%03d/group_ids_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	break;
  case 3:
	sprintf(buf, "%s/subhalo_ids_%03d",GrpPath,(int)Nsnap);
	break;
  case 4:
	sprintf(buf, "%s/group_ids_%03d",GrpPath,(int)Nsnap);
  break;
  default:
	fprintf(logfile,"error: grpfile_type not assigned? %s\n",buf);
	exit(1);
  }
  ByteOrder=check_grpcat_byteorder(buf);
  myfopen(fd,buf,"r");

  myfread(&Ngroups, sizeof(int), 1, fd);
  myfread(&TotNgroups, sizeof(int), 1, fd);
  Cat->Ngroups=TotNgroups;
  myfread(&Nids, sizeof(int), 1, fd);
  myfread(&TotNids,sizeof(long long),1,fd);
  #ifndef HBT_INT8
  if(TotNids>INT_MAX)
  {
	  printf("error: TotNids larger than HBTInt can hold! %lld\n",TotNids);
	  exit(1);
  }
  #endif
  Cat->Nids=TotNids;
  myfread(&NFiles, sizeof(int), 1, fd);
  //file offset:
  int dummy;
  myfread(&dummy,sizeof(int),1,fd);

  
  myfread(ICat.PID+Nload, sizeof(IDatInt), Nids, fd);
  if(feof(fd))
  {
	fprintf(logfile,"error:End-of-File in %s\n",buf);
	fflush(logfile);exit(1);  
  }
  Nload+=Nids;
  fclose(fd);
	}
	
  if(Nload!=Cat->Nids)
  {
	  fprintf(logfile,"error:Num grpparticles loaded not match: %lld,%lld\n",Nload,(long long)Cat->Nids);
	  fflush(logfile);
	  exit(1);
  }
  
  
#ifndef HBT_INT8
	  Cat->Len=ICat.Len;
	  Cat->Offset=ICat.Offset;
#else
	  Cat->Len=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  for(i=0;i<Cat->Ngroups;i++)
	  {
		  Cat->Len[i]=ICat.Len[i];
		  Cat->Offset[i]=ICat.Offset[i];
	  }
	  myfree(ICat.Len);
	  myfree(ICat.Offset);
#endif
  
  #ifdef HBTPID_RANKSTYLE
  IDatInt *PIDs,*p;  
  PIDs=load_PIDs_Sorted();
  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
  for(i=0;i<Cat->Nids;i++)
  {	  
	p=bsearch(&(ICat.PID[i]),PIDs,NP_DM,sizeof(IDatInt),comp_IDatInt);
	Cat->PIDorIndex[i]=p-PIDs;
  }
  myfree(PIDs);
  myfree(ICat.PID);
  #else
#ifdef SAME_INTTYPE
	  Cat->PIDorIndex=ICat.PID;
#else
	  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
	  for(i=0;i<Cat->Nids;i++)
	  Cat->PIDorIndex[i]=ICat.PID[i];
	  myfree(ICat.PID);
#endif 
  #endif
  
	for(i=1;i<Cat->Ngroups;i++)
	{
		if(Cat->Len[i]>Cat->Len[i-1])
		{
		printf("warning: groups not sorted with mass\n");
		}
	}
	
}	
void load_group_catalogue_v2(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{//old groupcat format
  FILE *fd;
  int i;
  char buf[1024];
  int TotNgroups,NFiles,dummy,Ngroups,Nids;
  int ByteOrder;
  
struct 
{
int *Len;
int *Offset;
IDatInt *PID; 
}ICat;

#ifdef SNAPLIST
Nsnap=snaplist[Nsnap];
#endif
  sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
  ByteOrder=check_grpcat_byteorder(buf);
  myfopen(fd,buf,"r");

  myfread(&Ngroups, sizeof(int), 1, fd);
  myfread(&Nids, sizeof(int), 1, fd);
  myfread(&TotNgroups, sizeof(int), 1, fd);
  myfread(&NFiles, sizeof(int), 1, fd);
  Cat->Ngroups=Ngroups;
  Cat->Nids=Nids;

  fprintf(logfile,"Snap="HBTIFMT"  Ngroups="HBTIFMT"  Nids="HBTIFMT"  TotNgroups=%d  NFiles=%d\n", Nsnap, Cat->Ngroups, Cat->Nids, TotNgroups, NFiles);

  ICat.Len= mymalloc(sizeof(int)*Cat->Ngroups);
  ICat.Offset=mymalloc(sizeof(int)*Cat->Ngroups);
  ICat.PID=mymalloc(sizeof(IDatInt)*Cat->Nids);
  Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);//consider move this out.............

  myfread(ICat.Len, sizeof(int), Cat->Ngroups, fd);
  myfread(ICat.Offset,sizeof(int), Cat->Ngroups, fd);
  fclose(fd);


  sprintf(buf, "%s/group_ids_%03d", GrpPath, (int)Nsnap);
  ByteOrder=check_grpcat_byteorder(buf);
  myfopen(fd,buf,"r");

  myfread(&dummy, sizeof(int), 1, fd);
  myfread(&dummy, sizeof(int), 1, fd);
  myfread(&dummy, sizeof(int), 1, fd);
  myfread(&dummy, sizeof(int), 1, fd);

  myfread(ICat.PID, sizeof(IDatInt), Cat->Nids, fd);
  
#ifndef HBT_INT8
	  Cat->Len=ICat.Len;
	  Cat->Offset=ICat.Offset;
#else
	  Cat->Len=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  for(i=0;i<Cat->Ngroups;i++)
	  {
		  Cat->Len[i]=ICat.Len[i];
		  Cat->Offset[i]=ICat.Offset[i];
	  }
	  myfree(ICat.Len);
	  myfree(ICat.Offset);
#endif
  
  #ifdef HBTPID_RANKSTYLE
  IDatInt *PIDs,*p;  
  PIDs=load_PIDs_Sorted();  
  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
  for(i=0;i<Cat->Nids;i++)
  {	  
	p=bsearch(&(ICat.PID[i]),PIDs,NP_DM,sizeof(IDatInt),comp_IDatInt);
	Cat->PIDorIndex[i]=p-PIDs;
  }
  myfree(PIDs);
  myfree(ICat.PID);
  #else
#ifdef SAME_INTTYPE
	  Cat->PIDorIndex=ICat.PID;
#else
	  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
	  for(i=0;i<Cat->Nids;i++)
	  Cat->PIDorIndex[i]=ICat.PID[i];
	  myfree(ICat.PID);
#endif
  #endif
  
  fclose(fd);
}
void load_group_catalogue_millimill(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{//milli-mellennium format
  FILE *fd;
  char buf[1024];
  int Ngroups,TotNgroups,Nids,NFiles;
  long long i,Nload,TotNids;
  int ByteOrder,grpfile_type;
  
struct 
{
int *Len;
int *Offset;
IDatInt *PID; 
}ICat;
  
  	#ifdef SNAPLIST
	Nsnap=snaplist[Nsnap];
	#endif

  Nload=0; TotNids=0;
  for(i=0;i<NFILES_GRP;i++)
  {
  grpfile_type=1;	
  sprintf(buf, "%s/groups_%03d/subhalo_tab_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
  if(!try_readfile(buf))
  {
	  sprintf(buf, "%s/snapdir_%03d/group_tab_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	  grpfile_type=2;
  }
  if(1==NFILES_GRP)
  {
	if(!try_readfile(buf))
	{
		sprintf(buf, "%s/subhalo_tab_%03d",GrpPath,(int)Nsnap);
		grpfile_type=3;
	}
	if(!try_readfile(buf))
	{
	sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
	grpfile_type=4;
	}
  }
	
  ByteOrder=check_grpcat_byteorder(buf);
  	
  myfopen(fd,buf,"r");

  myfread(&Ngroups, sizeof(int), 1, fd);
  myfread(&Nids, sizeof(int), 1, fd);
  myfread(&TotNgroups, sizeof(int), 1, fd);
  Cat->Ngroups=TotNgroups;
  TotNids+=Nids;
  //~ myfread(&TotNids,sizeof(long long),1,fd);
  myfread(&NFiles, sizeof(int), 1, fd);
 
  if(1==grpfile_type||3==grpfile_type)
  {
  int Nsub,TotNsub;
  myfread(&Nsub,sizeof(int),1,fd);
  myfread(&TotNsub,sizeof(int),1,fd);
  }
  
  if(NFILES_GRP!=NFiles)
	  {
		  fprintf(logfile,"error: number of grpfiles specified not the same as stored: %d,%d\n for file %s\n",
		  NFILES_GRP,(int)NFiles,buf);
		  fflush(logfile);
		  exit(1);
	  }
 
  if(0==i)
  fprintf(logfile,"Snap="HBTIFMT"  Ngroups=%d  Nids=%d  TotNgroups=%d  NFiles=%d\n", Nsnap, Ngroups, Nids, (int)(Cat->Ngroups), NFiles);
  
  if(0==i)
  {
  ICat.Len= mymalloc(sizeof(int)*Cat->Ngroups);
  ICat.Offset=mymalloc(sizeof(int)*Cat->Ngroups);
  }
  
  myfread(ICat.Len+Nload, sizeof(int), Ngroups, fd);
  myfread(ICat.Offset+Nload, sizeof(int), Ngroups, fd);
  if(feof(fd))
  {
	fprintf(logfile,"error:End-of-File in %s\n",buf);
	fflush(logfile);exit(1);  
  }
  Nload+=Ngroups;
  fclose(fd);
  }
  #ifndef HBT_INT8
  if(TotNids>INT_MAX)
  {
	  printf("error: TotNids larger than HBTInt can hold! %lld\n",TotNids);
	  exit(1);
  }
  #endif
  Cat->Nids=TotNids;
  if(Nload!=Cat->Ngroups)
  {
	  fprintf(logfile,"error:Num groups loaded not match: %lld,"HBTIFMT"\n",Nload,Cat->Ngroups);
	  fflush(logfile);exit(1);
  }
  
  ICat.PID=mymalloc(sizeof(IDatInt)*Cat->Nids);
  Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);//consider move this out.............
  Nload=0;
  for(i=0;i<NFILES_GRP;i++)
  {
  switch(grpfile_type)
  {
  case 1:
	sprintf(buf, "%s/groups_%03d/subhalo_ids_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	break;
  case 2:
	sprintf(buf, "%s/snapdir_%03d/group_ids_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	break;
  case 3:
	sprintf(buf, "%s/subhalo_ids_%03d",GrpPath,(int)Nsnap);
	break;
  case 4:
	sprintf(buf, "%s/group_ids_%03d",GrpPath,(int)Nsnap);
  break;
  default:
	fprintf(logfile,"error: grpfile_type not assigned? %s\n",buf);
	exit(1);
  }
  ByteOrder=check_grpcat_byteorder(buf);
  myfopen(fd,buf,"r");

  myfread(&Ngroups, sizeof(int), 1, fd);
  myfread(&Nids, sizeof(int), 1, fd);
  myfread(&TotNgroups, sizeof(int), 1, fd);
  myfread(&NFiles, sizeof(int), 1, fd);

  myfread(ICat.PID+Nload, sizeof(IDatInt), Nids, fd);
  if(feof(fd))
  {
	fprintf(logfile,"error:End-of-File in %s\n",buf);
	fflush(logfile);exit(1);  
  }
  Nload+=Nids;
  fclose(fd);
	}
	
  if(Nload!=Cat->Nids)
  {
	  fprintf(logfile,"error:Num grpparticles loaded not match: %lld,%lld\n",Nload,(long long)Cat->Nids);
	  fflush(logfile);
	  exit(1);
  }
  
  #define ID_MASK 0x00000003FFFFFFFF
  //~ if(Cat->Nids>0)
  //~ {
  //~ printf("%0llx,%0llx\n%lld,%lld\n",ICat.PID[0],ICat.PID[0]&ID_MASK,ICat.PID[0],ICat.PID[0]&ID_MASK);
  //~ exit(1);
  //~ }
  for(i=0;i<Cat->Nids;i++)  //only the least significant 34bits are the actual ids
	  ICat.PID[i]=ICat.PID[i]&ID_MASK;
  
#ifndef HBT_INT8
	  Cat->Len=ICat.Len;
	  Cat->Offset=ICat.Offset;
#else
	  Cat->Len=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  for(i=0;i<Cat->Ngroups;i++)
	  {
		  Cat->Len[i]=ICat.Len[i];
		  Cat->Offset[i]=ICat.Offset[i];
	  }
	  myfree(ICat.Len);
	  myfree(ICat.Offset);
#endif

TotNids=0;
	for(i=0;i<Cat->Ngroups;i++)
	{
		Cat->Offset[i]=TotNids;
		TotNids+=Cat->Len[i];
	}
	if(TotNids!=Cat->Nids)
   {
	  fprintf(logfile,"error:Num grpparticles summed not match: %lld,%lld\n",TotNids,(long long)Cat->Nids);
	  fflush(logfile);
	  exit(1);
   }
  
  #ifdef HBTPID_RANKSTYLE
  IDatInt *PIDs,*p;  
  PIDs=load_PIDs_Sorted();
  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
  for(i=0;i<Cat->Nids;i++)
  {	  
	p=bsearch(&(ICat.PID[i]),PIDs,NP_DM,sizeof(IDatInt),comp_IDatInt);
	Cat->PIDorIndex[i]=p-PIDs;
  }
  myfree(PIDs);
  myfree(ICat.PID);
  #else
#ifdef SAME_INTTYPE
	  Cat->PIDorIndex=ICat.PID;
#else
	  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
	  for(i=0;i<Cat->Nids;i++)
	  Cat->PIDorIndex[i]=ICat.PID[i];
	  myfree(ICat.PID);
#endif 
  #endif
  
	for(i=1;i<Cat->Ngroups;i++)
	{
		if(Cat->Len[i]>Cat->Len[i-1])
		{
		fprintf(logfile,"warning: groups not sorted with mass\n");
		}
	}
	
}	
void free_catalogue(CATALOGUE *A)
{
	free(A->Len);
	free(A->Offset);
	free(A->PIDorIndex);
	//~ free(A->HaloMask);
	//~ free(A->HaloMaskSrc);
	free(A->ID2Halo);
}
