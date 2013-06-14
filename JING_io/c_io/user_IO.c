#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

extern void open_fortran_file_(char *filename,int *fileno);
extern void read_fortran_record_(void *arr,int *nread,int *fileno);
extern void close_fortran_file_(int *fileno);

#ifdef BIGENDIAN
#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((unsigned) (x) >> 24) )
#define SWAP_8(x) ( ((x)<<56)|(((x) << 40)&0x00ff000000000000)|(((x)<<24)&0x0000ff0000000000)|(((x) << 8)&0x000000ff00000000)\
                 |(((x)>>8)&0x00000000ff000000)|(((x)>>24)&0x0000000000ff0000)|(((x)>>40)&0x000000000000ff00)|(((x)>>56)))
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))
#define FIX_DOUBLE(x) (*(unsigned long *)&(x) = SWAP_8(*(unsigned *)&(x)))
void swap_Nbyte(char *data,size_t n,size_t m)
/*This function is used to switch endian*/
{
  size_t i,j;
  char old_data[16];//by definition, sizeof(char)=1, one byte
  
  switch(m)
  {
	  case 1 :break;
	  case 2 :
	  		for(j=0;j<n;j++)
			FIX_SHORT(data[j*2]);
			break;
	  case 4 :
	  		for(j=0;j<n;j++)
			FIX_LONG(data[j*4]);
			break;
	  case 8 :
			for(j=0;j<n;j++)
			FIX_DOUBLE(data[j*8]);
			break;
	  default :
			for(j=0;j<n;j++)
			{
			  memcpy(&old_data[0],&data[j*m],m);
			  for(i=0;i<m;i++)
				{
				  data[j*m+i]=old_data[m-i-1];
				}
			}
	}
}
size_t fread_BE(void *buf,size_t Nsize,size_t Nbuf,FILE *fp)
{
	size_t Nread;
	Nread=fread(buf,Nsize,Nbuf,fp);
	swap_Nbyte((char *)buf,Nbuf,Nsize);
	return Nread;
}
#else
#define fread_BE fread
#endif
#ifndef PID_ORDERED
static int read_part_id_JING(int Nsnap,char *snapdir)
{
	FILE *fp;
	char buf[512];
	int i,dummy,dummy2;
	//~ size_t nread;
	#ifdef SKIP
		#undef SKIP
		#undef SKIP2
		#undef CHECK
	#endif
	#define SKIP fread_BE(&dummy,sizeof(dummy),1,fp)
	#define SKIP2 fread_BE(&dummy2,sizeof(dummy2),1,fp)
	#define CHECK if(dummy!=dummy2){fprintf(logfile,"error!record brackets not match for id%d!\t%d,%d\n",Nsnap,dummy,dummy2);fflush(logfile);exit(1);} 
	sprintf(buf,"%s/id%d.%04d",snapdir,RUN_NUM,Nsnap);
	if(!(fp=fopen(buf,"r")))
	{
		fprintf(logfile,"Error:cannot open file %s\n",buf);
		exit(1);
	}
	SKIP;
	fread_BE(Pdat.PID,sizeof(int),NP_DM,fp);
	SKIP2;
	CHECK;
	for(i=0;i<NP_DM;i++)
		if(Pdat.PID[i]<1||Pdat.PID[i]>NP_SIM)
			fprintf(logfile,"error: id not in the range 0~%d, for i=%d, pid=%d, snap=%d\n",NP_SIM,i,Pdat.PID[i],Nsnap);
	
	return NP_DM;
}
#endif
static int read_part_pos_JING(int Nsnap,char *snapdir)
{
	FILE *fp;
	char buf[512];
	short *tmp;
	int i,j,dummy,dummy2;
	//~ size_t nread;
	#ifdef SKIP
		#undef SKIP
		#undef SKIP2
		#undef CHECK
	#endif
	#define SKIP fread_BE(&dummy,sizeof(dummy),1,fp)
	#define SKIP2 fread_BE(&dummy2,sizeof(dummy2),1,fp)
	#define CHECK if(dummy!=dummy2){fprintf(logfile,"error!record brackets not match for pos%d!\t%d,%d\n",Nsnap,dummy,dummy2);fflush(logfile);exit(1);} 
	sprintf(buf,"%s/pos%d.%04d",snapdir,RUN_NUM,Nsnap);
	if(!(fp=fopen(buf,"r")))
	{
		fprintf(logfile,"Error:cannot open file %s\n",buf);
		exit(1);
	}
	#pragma omp critical(fill_header)
	{
	SKIP;
	fread_BE(&header.Np,sizeof(int),1,fp);
	fread_BE(&header.ips,sizeof(int),1,fp);
	fread_BE(&header.ztp,sizeof(float),1,fp);
	fread_BE(&header.Omegat,sizeof(float),1,fp);
	fread_BE(&header.Lambdat,sizeof(float),1,fp);
	fread_BE(&header.rLbox,sizeof(float),1,fp);
	fread_BE(&header.xscale,sizeof(float),1,fp);
	fread_BE(&header.vscale,sizeof(float),1,fp);
	SKIP2;
	CHECK;

	//~ printf("Np=%d,ips=%d,Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		//~ header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
	if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match %d\n",NP_DM);exit(1);}
	//~ 
	//~ printf("loading positions...\n");
	}
	if(Nsnap<=SNAP_DIV_SCALE)//has scale
	{
		tmp=malloc(sizeof(short)*header.Np);
		for(i=0;i<3;i++)
		{
			SKIP;
			fread_BE(tmp,sizeof(short),header.Np,fp);
			SKIP2;
			CHECK;
			for(j=0;j<header.Np;j++)
				Pdat.Pos[j][i]=tmp[j]*header.xscale;
		}
		free(tmp);
	}
	else
	{
		SKIP;
		for(i=0;i<3;i++)
			for(j=0;j<header.Np;j++)
	    	fread_BE(&Pdat.Pos[j][i],sizeof(float),1,fp);
		SKIP2;
		CHECK;
	}
	fclose(fp);
	return header.Np;
	#undef SKIP
	#undef SKIP2
	#undef CHECK
}
static int read_part_vel_JING(int Nsnap,char *snapdir)
{
	FILE *fp;
	char buf[512];
	short *tmp;
	int i,j,dummy,dummy2;
	#ifdef SKIP
		#undef SKIP
		#undef SKIP2
		#undef CHECK
	#endif
	#define SKIP fread_BE(&dummy,sizeof(dummy),1,fp)
	#define SKIP2 fread_BE(&dummy2,sizeof(dummy2),1,fp)
	#define CHECK if(dummy!=dummy2){fprintf(logfile,"error!record brackets not match for vel%d!\t%d,%d\n",Nsnap,dummy,dummy2);fflush(logfile);exit(1);} 

	sprintf(buf,"%s/vel%d.%04d",snapdir,RUN_NUM,Nsnap);
	if(!(fp=fopen(buf,"r")))
	{
		fprintf(logfile,"Error:cannot open file %s\n",buf);
		exit(1);
	}
	#pragma omp critical(fill_header)
	{
	SKIP;
	fread_BE(&header.Np,sizeof(int),1,fp);
	fread_BE(&header.ips,sizeof(int),1,fp);
	fread_BE(&header.ztp,sizeof(float),1,fp);
	fread_BE(&header.Omegat,sizeof(float),1,fp);
	fread_BE(&header.Lambdat,sizeof(float),1,fp);
	fread_BE(&header.rLbox,sizeof(float),1,fp);
	fread_BE(&header.xscale,sizeof(float),1,fp);
	fread_BE(&header.vscale,sizeof(float),1,fp);
	SKIP2;
	CHECK;
	//~ printf("Np=%d,ips=%d,Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		//~ header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
   if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match %d\n",NP_DM);exit(1);}
	//~ 
	//~ printf("loading velocities...\n");
	}
	if(Nsnap<=SNAP_DIV_SCALE)//has scale
	{
		tmp=malloc(sizeof(short)*header.Np);
		for(i=0;i<3;i++)
		{
			SKIP;
			fread_BE(tmp,sizeof(short),header.Np,fp);
			SKIP2;
			CHECK;
			for(j=0;j<header.Np;j++)
				Pdat.Vel[j][i]=tmp[j]*header.vscale;
		}
		free(tmp);
	}
	else
	{
		SKIP;
		for(i=0;i<3;i++)
			for(j=0;j<header.Np;j++)
			fread_BE(&Pdat.Vel[j][i],sizeof(float),1,fp);
		SKIP2;
		CHECK;
	}
	fclose(fp);
	return header.Np;
	#undef SKIP
	#undef SKIP2
	#undef CHECK
}	
void load_particle_data(int Nsnap,char *SnapPath)
{
	int i,j;
	float Hratio; //(Hz/H0)
	float scale_reduced,scale0;//a,R0

#pragma omp parallel sections
{
	#pragma omp section
	read_part_pos_JING(snaplist[Nsnap],SnapPath);
	#pragma omp section
	read_part_vel_JING(snaplist[Nsnap],SnapPath);
	#ifndef PID_ORDERED
	#pragma omp section
	read_part_id_JING(snaplist[Nsnap],SnapPath);
	#endif
}
	fprintf(logfile,"z=%g\n",header.ztp);fflush(logfile);
	Hratio=sqrt(OMEGAL0/header.Lambdat);
	header.Hz=HUBBLE0*Hratio;
	scale_reduced=1./(1.+header.ztp);
	header.time=scale_reduced;
	header.mass[0]=0.;
	header.mass[1]=PMass;
	scale0=1+Redshift_INI;//scale_INI=1,scale_reduced_INI=1./(1.+z_ini),so scale0=scale_INI/scale_reduce_INI;
	header.vunit=100.*header.rLbox*Hratio*scale_reduced*scale_reduced*scale0;   /*vunit=100*L*R*(H*R)/(H0*R0)
																				*      =100*L*(H/H0)*a*a*R0
																				* where a=R/R0;         */
	#pragma omp parallel for private(i,j)
	for(i=0;i<header.Np;i++)
		for(j=0;j<3;j++)
		{
			Pdat.Pos[i][j]-=floor(Pdat.Pos[i][j]);	//format coordinates to be in the range [0,1)
			Pdat.Pos[i][j]*=BOXSIZE;				//comoving coordinate in units of kpc/h
			Pdat.Vel[i][j]*=header.vunit;			//physical peculiar velocity in units of km/s
		}
}

void load_group_catalogue(int Nsnap,CATALOGUE *Cat,char *GrpPath)
{
  FILE *fp;
  char buf[1024];
  int i,dummy,dummy2;
  float b;
   #ifdef SKIP
		#undef SKIP
		#undef SKIP2
		#undef CHECK
	#endif
	#define SKIP fread_BE(&dummy,sizeof(dummy),1,fp)
	#define SKIP2 fread_BE(&dummy2,sizeof(dummy2),1,fp)
	#define CHECK if(dummy!=dummy2){fprintf(logfile,"error!record brackets not match for grp%d!\t%d,%d\n",Nsnap,dummy,dummy2);fflush(logfile);exit(1);} 

  sprintf(buf, "%s/fof.b20.%d.%04d",GrpPath,RUN_NUM,snaplist[Nsnap]);
  if(!(fp = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }
	
  SKIP;
  fread_BE(&b, sizeof(float),1,fp); // b is the link length
  fread_BE(&Cat->Ngroups, sizeof(int), 1, fp); //Nhalo is the number of halos in this snapshot
  SKIP2;
  CHECK;
  Cat->Nids=0;
  
  Cat->Len= mymalloc(sizeof(int)*Cat->Ngroups);
  Cat->Offset=mymalloc(sizeof(int)*Cat->Ngroups);
  Cat->HaloCen[0]=mymalloc(sizeof(float)*Cat->Ngroups);
  Cat->HaloCen[1]=mymalloc(sizeof(float)*Cat->Ngroups);
  Cat->HaloCen[2]=mymalloc(sizeof(float)*Cat->Ngroups);
  //~ Cat->HaloMask=mymalloc(sizeof(short)*NP_DM);
  Cat->ID2Halo=mymalloc(sizeof(int)*NP_DM);
  Cat->PIDorIndex=mymalloc(sizeof(int)*NP_DM);
  
    for(i=0;i<3;i++)
  {
    SKIP;
    fread_BE(Cat->HaloCen[i],sizeof(float),Cat->Ngroups,fp);
    SKIP2;
	CHECK;
  }
for(i=0;i<Cat->Ngroups;i++)
  {  
     
    SKIP;
    fread_BE(&Cat->Len[i],sizeof(int),1,fp);
    SKIP2;
	CHECK;
    Cat->HaloCen[0][i]*=BOXSIZE;
    Cat->HaloCen[1][i]*=BOXSIZE;
    Cat->HaloCen[2][i]*=BOXSIZE;
	Cat->Offset[i]=Cat->Nids;
	Cat->Nids+=Cat->Len[i];
	
    if(i>0)
    {
    	if(Cat->Len[i]>Cat->Len[i-1])
    	{
    	  fprintf(logfile,"wrong! in read_fof, %d, %d, %d\n",i,Cat->Len[i],Cat->Len[i-1]);
    	  exit(1);
    	}
    }
    
    SKIP;
    fread_BE(Cat->PIDorIndex+Cat->Offset[i],sizeof(int),Cat->Len[i],fp);
    SKIP2;
	CHECK;
  }
 
  Cat->PIDorIndex=realloc(Cat->PIDorIndex,sizeof(int)*Cat->Nids);
  fprintf(logfile,"Snap=%d SnapNum=%d  Ngroups=%d  Nids=%d b=%f\n",Nsnap,snaplist[Nsnap],Cat->Ngroups,Cat->Nids,b);

	for(i=0;i<Cat->Nids;i++)
	{
 	  if(Cat->PIDorIndex[i] == 0)
		fprintf(logfile,"i=%d groupPID=%d\n", i, (int)Cat->PIDorIndex[i]);//check if PID begin with ID=1;
	#ifdef GRPINPUT_INDEX	
		Cat->PIDorIndex[i]--;//change from [1,NP] to [0,NP-1] for index in C
	#endif
  	}

  fclose(fp);
  	#undef SKIP
	#undef SKIP2
	#undef CHECK
}
void free_catalogue(CATALOGUE *A)
{
	free(A->Len);
	free(A->HaloCen[0]);
	free(A->HaloCen[1]);
	free(A->HaloCen[2]);	
	free(A->Offset);
	free(A->PIDorIndex);
	//~ free(A->HaloMask);
	//~ free(A->HaloMaskSrc);
	free(A->ID2Halo);
}

