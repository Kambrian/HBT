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
extern void read_part_arr_imajor_(float partarr[][3],long int *np,int *fileno);//old format
extern void read_part_arr_xmajor_(float partarr[][3],long int *np,int *fileno);//newest format
extern void read_part_header_int4_(int *np,int *ips,float *ztp,float *omgt,float *lbdt,
			      float *boxsize,float *xscale,float *vscale,int *fileno);
extern void read_part_header_int8_(long int *np,long int *ips,float *ztp,float *omgt,float *lbdt,
			      float *boxsize,float *xscale,float *vscale,int *fileno);
extern void read_group_header_int4_(float *b,int *ngrp, int *fileno);
extern void read_group_header_int8_(float *b,long int *ngrp, int *fileno);

//Caution:: only support SAME_INTTYPE and SAME_REALTYPE; if they differ, need further modification of code!
#ifdef HBT_INT8
#define read_fortran_record_HBTInt read_fortran_record8_
#define read_part_header read_part_header_int8_
#define read_group_header read_group_header_int8_
#else
#define read_fortran_record_HBTInt read_fortran_record4_
#define read_part_header read_part_header_int4_
#define read_group_header read_group_header_int4_
#endif

#ifdef PDAT_XMAJOR
#define read_part_arr read_part_arr_xmajor_
#else
#define read_part_arr read_part_arr_imajor_
#endif

#ifdef BIGENDIAN
#define FLG_ENDIAN 1
#else
#define FLG_ENDIAN 0
#endif

#ifndef PID_ORDERED
void read_id_file_single(HBTInt *pid, long int np, int fileno, char *filename, int flag_endian)
{
	  int filestat;
	  open_fortran_file_(filename,&fileno,&flag_endian,&filestat);
	  if(filestat)
	  {
		  fprintf(logfile,"Error opening file %s,error no. %d\n",filename,filestat);fflush(logfile);
		  exit(1);
	  }
	  read_fortran_record_HBTInt(pid, &np, &fileno);
	  close_fortran_file_(&fileno);
}
static HBTInt read_part_id_JING(int Nsnap,char *snapdir)
{
	char buf[1024];
	HBTInt i;
	long int nread;
	int fileno=10, ifile;
	
	Pdat.PID=mymalloc(sizeof(HBTInt)*NP_DM);
		
#if !defined(NFILE_ID)||NFILE_ID==1
	sprintf(buf,"%s/id%d.%04d",snapdir,RUN_NUM,(int)Nsnap);
	read_id_file_single(Pdat.PID, NP_DM, fileno, buf, FLG_ENDIAN);
#else
	nread=NP_DM/NFILE_ID;
	#pragma omp parallel for num_threads(NFILE_ID) private(ifile,buf)
	for(ifile=0;ifile<NFILE_ID;ifile++)
	{
	  sprintf(buf,"%s/id%d.%04d.%02d",snapdir,RUN_NUM,(int)Nsnap, ifile+1);
	  read_id_file_single(Pdat.PID+nread*ifile, nread, fileno+ifile, buf, FLG_ENDIAN);
	}
#endif	
#pragma omp parallel for
	for(i=0;i<NP_DM;i++)
		if(Pdat.PID[i]<1||Pdat.PID[i]>NP_SIM)
		{fprintf(logfile,"error: id not in the range 0~"HBTIFMT", for i="HBTIFMT", pid="HBTIFMT", snap=%d\n",NP_SIM,i,Pdat.PID[i],Nsnap);fflush(logfile);}
	
	return NP_DM;
}
#endif
void read_pos_file_single(HBTReal pos[][3], long int np, int fileno, char *filename, int flag_endian, int has_header, int has_scale)
{
  int filestat;
  HBTInt i,j;
  open_fortran_file_(filename,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",filename,filestat);fflush(logfile);
		exit(1);
	}
	
	if(has_header)
	{
	read_part_header(&header.Np,&header.ips,&header.ztp,&header.Omegat,&header.Lambdat,
						&header.rLbox,&header.xscale,&header.vscale,&fileno);
	 printf("Np="HBTIFMT",ips="HBTIFMT",Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		 header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
	if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match in %s\n "HBTIFMT"\n", filename, header.Np);fflush(logfile);exit(1);}
	}
	
// 	Pdat.Pos=mymalloc(sizeof(HBTReal)*3*NP_DM);
	if(has_scale)//has scale
	{
	  short *tmp;
		if(sizeof(short)!=2)
		{
			fprintf(logfile,"error: sizeof(short)=%zd, !=2\n",sizeof(short));fflush(logfile);
			exit(1);
		}
		tmp=malloc(sizeof(short)*np);
		for(i=0;i<3;i++)
		{
			read_fortran_record2_(tmp,&np,&fileno);
			for(j=0;j<np;j++)
				pos[j][i]=tmp[j]*header.xscale;
		}
		free(tmp);
	}
	else
		read_part_arr(pos,&np,&fileno);
	close_fortran_file_(&fileno);
}
void read_vel_file_single(HBTReal vel[][3], long int np, int fileno, char *filename, int flag_endian, int has_header, int has_scale)
{
  int filestat;
  HBTInt i,j;
  open_fortran_file_(filename,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",filename,filestat);fflush(logfile);
		exit(1);
	}
	
	if(has_header)
	{
	read_part_header(&header.Np,&header.ips,&header.ztp,&header.Omegat,&header.Lambdat,
						&header.rLbox,&header.xscale,&header.vscale,&fileno);
	 printf("Np="HBTIFMT",ips="HBTIFMT",Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		 header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
	if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match in %s\n "HBTIFMT"\n", filename, header.Np);fflush(logfile);exit(1);}
	}
	
	if(has_scale)//has scale
	{
	  short *tmp;
		if(sizeof(short)!=2)
		{
			fprintf(logfile,"error: sizeof(short)=%zd, !=2\n",sizeof(short));fflush(logfile);
			exit(1);
		}
		tmp=malloc(sizeof(short)*np);
		for(i=0;i<3;i++)
		{
			read_fortran_record2_(tmp,&np,&fileno);
			for(j=0;j<np;j++)
				vel[j][i]=tmp[j]*header.vscale;
		}
		free(tmp);
	}
	else
		read_part_arr(vel,&np,&fileno);
	close_fortran_file_(&fileno);
}
static HBTInt read_part_pos_JING(int Nsnap,char *snapdir)
{
	char buf[1024];
	int fileno=11,ifile;
	long int nread;
	
	Pdat.Pos=mymalloc(sizeof(HBTReal)*3*NP_DM);
	
#if !defined(NFILE_POS)||NFILE_POS==1
	sprintf(buf,"%s/pos%d.%04d",snapdir,RUN_NUM,Nsnap);
	read_pos_file_single(Pdat.Pos, NP_DM, fileno, buf, FLG_ENDIAN, 1, Nsnap<=SNAP_DIV_SCALE);
#else
	nread=NP_DM/NFILE_POS;
	#pragma omp parallel for num_threads(NFILE_POS) private(ifile,buf)
	for(ifile=0;ifile<NFILE_POS;ifile++)
	{
	  sprintf(buf,"%s/pos%d.%04d.%02d",snapdir,RUN_NUM,(int)Nsnap, ifile+1);
	  read_pos_file_single(Pdat.Pos+nread*ifile, nread, fileno+ifile, buf, FLG_ENDIAN, 0==ifile,  Nsnap<=SNAP_DIV_SCALE);
	}
#endif	
	return header.Np;
}
static HBTInt read_part_vel_JING(int Nsnap,char *snapdir)
{
	char buf[1024];
	int fileno=12,ifile;
	long int nread;
	
	Pdat.Vel=mymalloc(sizeof(HBTReal)*NP_DM*3);
#if !defined(NFILE_VEL)||NFILE_VEL==1	
	sprintf(buf,"%s/vel%d.%04d",snapdir,RUN_NUM,Nsnap);
	read_vel_file_single(Pdat.Vel, NP_DM, fileno, buf, FLG_ENDIAN, 1, Nsnap<=SNAP_DIV_SCALE);
#else
	nread=NP_DM/NFILE_VEL;
	#pragma omp parallel for num_threads(NFILE_VEL) private(ifile,buf)
	for(ifile=0;ifile<NFILE_VEL;ifile++)
	{
	  sprintf(buf,"%s/vel%d.%04d.%02d",snapdir,RUN_NUM,(int)Nsnap, ifile+1);
	  read_vel_file_single(Pdat.Vel+nread*ifile, nread, fileno+ifile, buf, FLG_ENDIAN, 0==ifile,  Nsnap<=SNAP_DIV_SCALE);
	}
#endif
	return header.Np;
}

void load_particle_header(HBTInt Nsnap, char *SnapPath)
{	
	float Hratio; //(Hz/H0)
	float scale_reduced,scale0;//a,R0

	char buf[1024];
	short *tmp;
	int flag_endian,filestat,fileno=13;
	
	flag_endian=FLG_ENDIAN;
#if !defined(NFILE_POS)||NFILE_POS==1	
	sprintf(buf,"%s/pos%d.%04d",SnapPath,RUN_NUM,(int)snaplist[Nsnap]);
#else
	sprintf(buf,"%s/pos%d.%04d.%02d",SnapPath,RUN_NUM,(int)snaplist[Nsnap], 1);
#endif
	open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",buf,filestat);fflush(logfile);
		exit(1);
	}

//	#pragma omp critical (fill_header)
	{
	read_part_header(&header.Np,&header.ips,&header.ztp,&header.Omegat,&header.Lambdat,
						&header.rLbox,&header.xscale,&header.vscale,&fileno);
	//~ printf("Np=%d,ips=%d,Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\n",
		//~ header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale);
	if(header.Np!=NP_DM){fprintf(logfile,"error:particle number do not match(pos) "HBTIFMT"\n",header.Np);fflush(logfile);exit(1);}
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
	header.vunit=100.*header.rLbox*Hratio*scale_reduced*scale_reduced*scale0;   /*vunit=100*rLbox*R*(H*R)/(H0*R0)
											  =L*H0*Hratio*R*R/R0 (H0=100 when length(L) in Mpc/h)
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
	HBTInt i,j;
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
void load_particle_data(HBTInt Nsnap,char *SnapPath)
{
	HBTInt i,j;
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
void load_group_catalogue(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{
	load_group_catalogue_HBT(Nsnap,Cat,GrpPath);	
}
#else
void load_group_catalogue(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{
  char buf[1024];
  float b;
  HBTInt i,j;
  int flag_endian,filestat,fileno=13;
  long int nread;
	
	flag_endian=FLG_ENDIAN;
    sprintf(buf, "%s/fof.b20.%d.%04d",GrpPath,RUN_NUM,(int)snaplist[Nsnap]);
	open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
	if(filestat)
	{
		fprintf(logfile,"Error opening file %s,error no. %d\n",buf,filestat);fflush(logfile);
		exit(1);
	}
	read_group_header(&b,&Cat->Ngroups,&fileno);	

  Cat->Nids=0;
  
  Cat->Len= mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  Cat->HaloCen[0]=mymalloc(sizeof(float)*Cat->Ngroups);
  Cat->HaloCen[1]=mymalloc(sizeof(float)*Cat->Ngroups);
  Cat->HaloCen[2]=mymalloc(sizeof(float)*Cat->Ngroups);
  //~ Cat->HaloMask=mymalloc(sizeof(short)*NP_DM);
  Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);
  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*NP_DM);
  
  nread=Cat->Ngroups;
  for(i=0;i<3;i++)
	read_fortran_record4_(Cat->HaloCen[i],&nread,&fileno);

for(i=0;i<Cat->Ngroups;i++)
  {  
    nread=1;
	read_fortran_record_HBTInt(&Cat->Len[i],&nread,&fileno); 
    Cat->HaloCen[0][i]*=BOXSIZE;
    Cat->HaloCen[1][i]*=BOXSIZE;
    Cat->HaloCen[2][i]*=BOXSIZE;
	Cat->Offset[i]=Cat->Nids;
	Cat->Nids+=Cat->Len[i];
	
    if(i>0)
    {
    	if(Cat->Len[i]>Cat->Len[i-1])
    	{
    	  fprintf(logfile,"wrong! in read_fof, "HBTIFMT", "HBTIFMT", "HBTIFMT"\n",i,Cat->Len[i],Cat->Len[i-1]);fflush(logfile);
    	  exit(1);
    	}
    }
	
    nread=Cat->Len[i];
	read_fortran_record_HBTInt(Cat->PIDorIndex+Cat->Offset[i],&nread,&fileno); 
  }
 
 close_fortran_file_(&fileno);
 
  Cat->PIDorIndex=realloc(Cat->PIDorIndex,sizeof(HBTInt)*Cat->Nids);
  fprintf(logfile,"Snap="HBTIFMT" SnapNum="HBTIFMT"  Ngroups="HBTIFMT"  Nids="HBTIFMT" b=%f\n",Nsnap,snaplist[Nsnap],Cat->Ngroups,Cat->Nids,b);

	for(i=0;i<Cat->Nids;i++)
	{
 	  if(Cat->PIDorIndex[i] == 0)
		fprintf(logfile,"i="HBTIFMT" groupPID="HBTIFMT"\n", i, Cat->PIDorIndex[i]);//check if PID begin with ID=1;
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