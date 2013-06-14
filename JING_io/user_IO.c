#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

extern void open_fortran_file_(char *filename,int *fileno,int *endian,int *error);
extern void read_fortran_record1_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record2_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record4_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record8_(void *arr,long int *arr_len,int *fileno);
extern void close_fortran_file_(int *fileno);
extern void read_part_arr_(float partarr[][3],int *np,int *fileno);
extern void read_part_header_(int *np,int *ips,float *ztp,float *omgt,float *lbdt,
									float *boxsize,float *xscale,float *vscale,int *fileno);
#ifdef BIGENDIAN
#define FLG_ENDIAN 1
#else
#define FLG_ENDIAN 0
#endif

#ifndef PID_ORDERED
static int read_part_id_JING(int Nsnap,char *snapdir)
{
	char buf[1024];
	int i,flag_endian,filestat,fileno=10;
	long int nread;
	
	flag_endian=FLG_ENDIAN;
	sprintf(buf,"%s/id%d.%04d",snapdir,RUN_NUM,Nsnap);
	open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",buf,filestat);fflush(logfile);
		exit(1);
	}
	nread=NP_DM;
	Pdat.PID=mymalloc(sizeof(int)*NP_DM);
	read_fortran_record4_(Pdat.PID,&nread,&fileno);	
	for(i=0;i<NP_DM;i++)
		if(Pdat.PID[i]<1||Pdat.PID[i]>NP_SIM)
			fprintf(logfile,"error: id not in the range 0~%d, for i=%d, pid=%d, snap=%d\n",NP_SIM,i,Pdat.PID[i],Nsnap);fflush(logfile);
	close_fortran_file_(&fileno);
	return NP_DM;
}
#endif
static int read_part_pos_JING(int Nsnap,char *snapdir)
{
	char buf[1024];
	short *tmp;
	int i,j,flag_endian,filestat,fileno=11;
	long int nread;
	
	flag_endian=FLG_ENDIAN;		
	sprintf(buf,"%s/pos%d.%04d",snapdir,RUN_NUM,Nsnap);
	open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",buf,filestat);fflush(logfile);
		exit(1);
	}

//	#pragma omp critical (fill_header)
	{
	read_part_header_(&header.Np,&header.ips,&header.ztp,&header.Omegat,&header.Lambdat,
						&header.rLbox,&header.xscale,&header.vscale,&fileno);
	//~ printf("Np=%d,ips=%d,Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		//~ header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
	if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match(pos) %d\n",NP_DM);fflush(logfile);exit(1);}
	}
	Pdat.Pos=mymalloc(sizeof(float)*3*NP_DM);
	if(Nsnap<=SNAP_DIV_SCALE)//has scale
	{
		nread=header.Np;
		if(sizeof(short)!=2)
		{
			fprintf(logfile,"error: sizeof(short)=%ld, !=2\n",sizeof(short));fflush(logfile);
			exit(1);
		}
		tmp=malloc(sizeof(short)*nread);
		for(i=0;i<3;i++)
		{
			read_fortran_record2_(tmp,&nread,&fileno);
			for(j=0;j<header.Np;j++)
				Pdat.Pos[j][i]=tmp[j]*header.xscale;
		}
		free(tmp);
	}
	else
	{
		//~ nread=3*header.Np;
		//~ read_fortran_record4_(&(Pdat.Pos[0][0]),&nread,&fileno);
		read_part_arr_(Pdat.Pos,&header.Np,&fileno);
	}
	close_fortran_file_(&fileno);
	return header.Np;
}
static int read_part_vel_JING(int Nsnap,char *snapdir)
{
	char buf[1024];
	short *tmp;
	int i,j,flag_endian,filestat,fileno=12;
	long int nread;
	
	flag_endian=FLG_ENDIAN;
	sprintf(buf,"%s/vel%d.%04d",snapdir,RUN_NUM,Nsnap);
	open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",buf,filestat);fflush(logfile);
		exit(1);
	}

//	#pragma omp critical (fill_header)
	{	
	read_part_header_(&header.Np,&header.ips,&header.ztp,&header.Omegat,&header.Lambdat,
						&header.rLbox,&header.xscale,&header.vscale,&fileno);
	//~ printf("Np=%d,ips=%d,Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		//~ header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
	if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match(vel) %d\n",NP_DM);fflush(logfile);exit(1);}
	}
	
	Pdat.Vel=mymalloc(sizeof(float)*3*NP_DM);
	if(Nsnap<=SNAP_DIV_SCALE)//has scale
	{
		nread=header.Np;
		if(sizeof(short)!=2)
		{
			fprintf(logfile,"error: sizeof(short)=%ld, !=2\n",sizeof(short));fflush(logfile);
			exit(1);
		}
		tmp=malloc(sizeof(short)*nread);
		for(i=0;i<3;i++)
		{
			read_fortran_record2_(tmp,&nread,&fileno);
			for(j=0;j<header.Np;j++)
				Pdat.Vel[j][i]=tmp[j]*header.vscale;
		}
		free(tmp);
	}
	else
	{
		//~ nread=3*header.Np;
		//~ read_fortran_record4_(&(Pdat.Vel[0][0]),&nread,&fileno);
		read_part_arr_(Pdat.Vel,&header.Np,&fileno);
	}
	close_fortran_file_(&fileno);
	return header.Np;
}

void load_particle_header(int Nsnap, char *SnapPath)
{	
	float Hratio; //(Hz/H0)
	float scale_reduced,scale0;//a,R0

	char buf[1024];
	short *tmp;
	int flag_endian,filestat,fileno=13;
	
	flag_endian=FLG_ENDIAN;		
	sprintf(buf,"%s/pos%d.%04d",SnapPath,RUN_NUM,snaplist[Nsnap]);
	open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",buf,filestat);fflush(logfile);
		exit(1);
	}

//	#pragma omp critical (fill_header)
	{
	read_part_header_(&header.Np,&header.ips,&header.ztp,&header.Omegat,&header.Lambdat,
						&header.rLbox,&header.xscale,&header.vscale,&fileno);
	//~ printf("Np=%d,ips=%d,Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		//~ header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
	if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match(pos) %d\n",NP_DM);fflush(logfile);exit(1);}
	}
	close_fortran_file_(&fileno);

	fprintf(logfile,"z=%g\n",header.ztp);fflush(logfile);
	Hratio=sqrt(OMEGAL0/header.Lambdat);
	header.Hz=HUBBLE0*Hratio;
	scale_reduced=1./(1.+header.ztp);
	header.time=scale_reduced;
	header.mass[0]=0.;
	header.mass[1]=PMass;
	header.Omega0=OMEGA0;
	header.OmegaLambda=OMEGAL0;	
	scale0=1+Redshift_INI;//scale_INI=1,scale_reduced_INI=1./(1.+z_ini),so scale0=scale_INI/scale_reduce_INI;
	header.vunit=100.*header.rLbox*Hratio*scale_reduced*scale_reduced*scale0;   /*vunit=100*L*R*(H*R)/(H0*R0)
																				*      =100*L*(H/H0)*a*a*R0
																				* where a=R/R0;         */
	header.Nsnap=Nsnap;
}


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files). use the bitwise loadflags to 
 * determine which part of the data to load. the last 3 bits
 * of loadflags encodes whether to load vel, pos and id.
 */
void load_particle_data_bypart(HBTInt Nsnap, char *SnapPath, unsigned char loadflags)
{
	int i,j;
	float Hratio; //(Hz/H0)
	float scale_reduced,scale0;//a,R0
	unsigned char flag_id,flag_pos,flag_vel;
	
	/*loadflags encodes the flag to load id,pos,vel in its lowest,second lowest and third lowest bit */
	flag_id=get_bit(loadflags,0);
	flag_pos=get_bit(loadflags,1);
	flag_vel=get_bit(loadflags,2);
	
	Pdat.Nsnap=Nsnap;
	
	if(flag_pos)
		read_part_pos_JING(snaplist[Nsnap],SnapPath);
	else
		Pdat.Pos=NULL;	
	
	if(flag_vel)
		read_part_vel_JING(snaplist[Nsnap],SnapPath);
	else
		Pdat.Vel=NULL;	
	
	if(flag_id)
	{
	#ifndef PID_ORDERED
		read_part_id_JING(snaplist[Nsnap],SnapPath);
	#endif
	}
	else
	{
		Pdat.PID=NULL;
	}
	
	fprintf(logfile,"z=%g\n",header.ztp);fflush(logfile);
	Hratio=sqrt(OMEGAL0/header.Lambdat);
	header.Hz=HUBBLE0*Hratio;
	scale_reduced=1./(1.+header.ztp);
	header.time=scale_reduced;
	header.mass[0]=0.;
	header.mass[1]=PMass;
	header.Omega0=OMEGA0;
	header.OmegaLambda=OMEGAL0;	
	scale0=1+Redshift_INI;//scale_INI=1,scale_reduced_INI=1./(1.+z_ini),so scale0=scale_INI/scale_reduce_INI;
	header.vunit=100.*header.rLbox*Hratio*scale_reduced*scale_reduced*scale0;   /*vunit=100*L*R*(H*R)/(H0*R0)
																				*      =100*L*(H/H0)*a*a*R0
																				* where a=R/R0;         */
	header.Nsnap=Nsnap;
	if(flag_pos)
		#pragma omp parallel for private(i,j)
		for(i=0;i<header.Np;i++)
			for(j=0;j<3;j++)
			{
				Pdat.Pos[i][j]-=floor(Pdat.Pos[i][j]);	//format coordinates to be in the range [0,1)
				Pdat.Pos[i][j]*=BOXSIZE;				//comoving coordinate in units of kpc/h
			}
	if(flag_vel)	
		#pragma omp parallel for private(i,j)
		for(i=0;i<header.Np;i++)
			for(j=0;j<3;j++)
				Pdat.Vel[i][j]*=header.vunit;			//physical peculiar velocity in units of km/s

}
void load_particle_data(int Nsnap,char *SnapPath)
{
	int i,j;
	float Hratio; //(Hz/H0)
	float scale_reduced,scale0;//a,R0
//	#pragma omp parallel sections
	{
//	#pragma omp section
	read_part_pos_JING(snaplist[Nsnap],SnapPath);
//	#pragma omp section
	read_part_vel_JING(snaplist[Nsnap],SnapPath);
	#ifndef PID_ORDERED
//	#pragma omp section
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
	header.Omega0=OMEGA0;
	header.OmegaLambda=OMEGAL0;	
	scale0=1+Redshift_INI;//scale_INI=1,scale_reduced_INI=1./(1.+z_ini),so scale0=scale_INI/scale_reduce_INI;
	header.vunit=100.*header.rLbox*Hratio*scale_reduced*scale_reduced*scale0;   /*vunit=100*L*R*(H*R)/(H0*R0)
																				*      =100*L*(H/H0)*a*a*R0
																				* where a=R/R0;         */
	header.Nsnap=Nsnap;
	#pragma omp parallel for private(i,j)
	for(i=0;i<header.Np;i++)
		for(j=0;j<3;j++)
		{
			Pdat.Pos[i][j]-=floor(Pdat.Pos[i][j]);	//format coordinates to be in the range [0,1)
			Pdat.Pos[i][j]*=BOXSIZE;				//comoving coordinate in units of kpc/h
			Pdat.Vel[i][j]*=header.vunit;			//physical peculiar velocity in units of km/s
		}
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
#ifdef GRP_HBTFORMAT
void load_group_catalogue(int Nsnap,CATALOGUE *Cat,char *GrpPath)
{
	load_group_catalogue_HBT(Nsnap,Cat,GrpPath);	
}
#else
void load_group_catalogue(int Nsnap,CATALOGUE *Cat,char *GrpPath)
{
  char buf[1024];
  float b,header_arr[2];
  int i,j,flag_endian,filestat,fileno=13;
  long int nread;
	
	flag_endian=FLG_ENDIAN;
    sprintf(buf, "%s/fof.b20.%d.%04d",GrpPath,RUN_NUM,snaplist[Nsnap]);
	open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",buf,filestat);fflush(logfile);
		exit(1);
	}
	nread=2;
	read_fortran_record4_(header_arr,&nread,&fileno);	
	b=header_arr[0];
	Cat->Ngroups=*((int *)(header_arr+1));

  Cat->Nids=0;
  
  Cat->Len= mymalloc(sizeof(int)*Cat->Ngroups);
  Cat->Offset=mymalloc(sizeof(int)*Cat->Ngroups);
  Cat->HaloCen[0]=mymalloc(sizeof(float)*Cat->Ngroups);
  Cat->HaloCen[1]=mymalloc(sizeof(float)*Cat->Ngroups);
  Cat->HaloCen[2]=mymalloc(sizeof(float)*Cat->Ngroups);
  //~ Cat->HaloMask=mymalloc(sizeof(short)*NP_DM);
  Cat->ID2Halo=mymalloc(sizeof(int)*NP_DM);
  Cat->PIDorIndex=mymalloc(sizeof(int)*NP_DM);
  
  nread=Cat->Ngroups;
  for(i=0;i<3;i++)
	read_fortran_record4_(Cat->HaloCen[i],&nread,&fileno);

for(i=0;i<Cat->Ngroups;i++)
  {  
    nread=1;
	read_fortran_record4_(&Cat->Len[i],&nread,&fileno); 
    Cat->HaloCen[0][i]*=BOXSIZE;
    Cat->HaloCen[1][i]*=BOXSIZE;
    Cat->HaloCen[2][i]*=BOXSIZE;
	Cat->Offset[i]=Cat->Nids;
	Cat->Nids+=Cat->Len[i];
	
    if(i>0)
    {
    	if(Cat->Len[i]>Cat->Len[i-1])
    	{
    	  fprintf(logfile,"wrong! in read_fof, %d, %d, %d\n",i,Cat->Len[i],Cat->Len[i-1]);fflush(logfile);
    	  exit(1);
    	}
    }
	
    nread=Cat->Len[i];
	read_fortran_record4_(Cat->PIDorIndex+Cat->Offset[i],&nread,&fileno); 
  }
 
 close_fortran_file_(&fileno);
 
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
}
#endif
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

#undef FLG_ENDIAN
