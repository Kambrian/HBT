#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

struct BigData
{
int *PID; 
float (*Pos)[3];
float (*Vel)[3];
float *Mass;
} Bdat;	

void load_big_particles(int Nsnap, char *SnapPath);
void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);

#define NGRID 1500
	int mapxy[NGRID][NGRID]={0},mapxz[NGRID][NGRID]={0},mapyz[NGRID][NGRID]={0};
	
int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat;
	int Nsnap,i,j,pid,subid,grpid,nmass;
	FILE *fp,*fpr,*fpv;
	char buf[1024];
	float step[3],range[3][2],rvir;
	float (*Mvir)[3],(*Rvir)[3],*cen,bmass;

	int grid[3];
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/image",SUBCAT_DIR);
	logfile=stdout;

	if(argc!=3)
	{
		printf("usage: %s [Nsnap] [grpid]\n",argv[0]);
		exit(1);
	}
	else
	{
		Nsnap=atoi(argv[1]);
		grpid=atoi(argv[2]);
	}

//==================fof 2d map======================//	
	load_big_particles(Nsnap,SNAPSHOT_DIR);
	//~ rvir=comoving_virial_radius(Cat.Len[grpid]);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,header.mass[1],SubCat.Ngroups,Nsnap);
	rvir=Rvir[grpid][0];
	subid=SubCat.GrpOffset_Sub[grpid];
	cen=SubCat.Property[subid].CoM;
	for(i=0;i<3;i++)
	{
			range[i][0]=cen[i]-rvir;
			range[i][1]=cen[i]+rvir;
	}
	sprintf(buf,"%s/fofbsize_%03d_%d",outputdir,Nsnap,grpid);
	myfopen(fp,buf,"w");
	for(i=0;i<3;i++)	
		fprintf(fp,"%f,%f\n",range[i][0],range[i][1]);
	fclose(fp);
	
	for(i=0;i<3;i++)
		step[i]=(range[i][1]-range[i][0])/NGRID;
	for(i=0;i<NGRID;i++)
		for(j=0;j<NGRID;j++)
		{
			mapxy[i][j]=0;
			mapxz[i][j]=0;
			mapyz[i][j]=0;
		}
	bmass=0;nmass=0;		
	for(i=0;i<NP_SIM-NP_DM-NP_gas;i++)
	{
		pid=i;
		if(distance(SubCat.Property[subid].CoM,Bdat.Pos[pid])>=rvir) continue;
		for(j=0;j<3;j++)
		{
			grid[j]=floor((Bdat.Pos[pid][j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
			if(grid[j]<0) grid[j]=0;
		}
		mapxy[grid[0]][grid[1]]++;
		mapxz[grid[0]][grid[2]]++;
		mapyz[grid[1]][grid[2]]++;
		nmass++;
		bmass+=Bdat.Mass[pid];
	}
	
	printf("Contamination fraction: %g (%d particles)\n",bmass/Mvir[grpid][0],nmass);
	
	sprintf(buf,"%s/fofbmapxy_%03d_%d",outputdir,Nsnap,grpid);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxy[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/fofbmapxz_%03d_%d",outputdir,Nsnap,grpid);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/fofbmapyz_%03d_%d",outputdir,Nsnap,grpid);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapyz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	
/*===============================box 2d map=======================================*/
	for(i=0;i<3;i++)
	{
			range[i][0]=0.;
			range[i][1]=BOXSIZE;
	}
	
	for(i=0;i<3;i++)
		step[i]=(range[i][1]-range[i][0])/NGRID;
	
	sprintf(buf,"%s/bigsize_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<3;i++)	
		fprintf(fp,"%f,%f\n",range[i][0],range[i][1]);
	fclose(fp);

	for(i=0;i<NGRID;i++)
		for(j=0;j<NGRID;j++)
		{
			mapxy[i][j]=0;
			mapxz[i][j]=0;
			mapyz[i][j]=0;
		}
	for(i=0;i<NP_SIM-NP_DM-NP_gas;i++)
	{
		pid=i;
		for(j=0;j<3;j++)
		{
			grid[j]=floor((Bdat.Pos[pid][j]-range[j][0])/step[j]);
			if(grid[j]>=NGRID) grid[j]=NGRID-1;
		}
		mapxy[grid[0]][grid[1]]++;
		mapxz[grid[0]][grid[2]]++;
		mapyz[grid[1]][grid[2]]++;
	}
	
	sprintf(buf,"%s/bigmapxy_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxy[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/bigmapxz_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapxz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	//~ 
	sprintf(buf,"%s/bigmapyz_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
			fprintf(fp,"%d\t",mapyz[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	myfree(Bdat.Mass);
	myfree(Bdat.Pos);
return 0;
}

		
void load_big_particles(int Nsnap, char *SnapPath)
{
	int Ndm,Ngas,Nother,Nmass,i;
	int dummy,dummy2;
	long int pre_len,tail_len;
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


	/*=====reading snapshot========*/
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
	//fprintf(logfile,"headerOK\n");
	header.Hz=HUBBLE0 * sqrt(header.Omega0 / (header.time * header.time * header.time) 
			+ (1 - header.Omega0 - header.OmegaLambda) / (header.time * header.time)
			+ header.OmegaLambda);//Hubble param for the current catalogue;

	Ngas=header.npart[0];	Ndm=header.npart[1];	Nother=header.npart[2]+header.npart[3]+header.npart[4]+header.npart[5];
	for(Nmass=0,i=0;i<6;i++)
	{
		if(!header.mass[i]) Nmass+=header.npart[i];
	}
	if(Nother!=Nmass)
	{
		fprintf(logfile,"error: Nother=%d,Nmass=%d\n",Nother,Nmass);
		exit(1);
	}
	//Bdat.PID=mymalloc(sizeof(int)*Nother);		
	Bdat.Pos=mymalloc(sizeof(float)*3*Nother);
	//Bdat.Vel=mymalloc(sizeof(float)*3*Nother);
	Bdat.Mass=mymalloc(sizeof(float)*Nother);
	
	pre_len=(NP_gas+NP_DM)*sizeof(float)*3;
	
	SKIP;
	fseek(fp,pre_len,SEEK_CUR);
	fread(Bdat.Pos,sizeof(float)*3,Nother,fp);
	SKIP2;
	CHECK;
	//fprintf(logfile,"PosOK\n");
	
	SKIP;
	fseek(fp,sizeof(float)*3*NP_SIM,SEEK_CUR);
	//~ fseek(fp,pre_len,SEEK_CUR);
	//~ fread(Bdat.Vel,sizeof(float)*3,Nother,fp);
	SKIP2;
	CHECK;
	//fprintf(logfile,"VelOK\n");
	
	
	SKIP;
	fseek(fp,sizeof(int)*NP_SIM,SEEK_CUR);
	//~ pre_len=(NP_gas+NP_DM)*sizeof(int);
	//~ fseek(fp,pre_len,SEEK_CUR);
	//~ fread(Bdat.PID,sizeof(int),Nother,fp);
	SKIP2;
	CHECK;
	//fprintf(logfile,"IDOK\n");
	
	SKIP;
	fread(Bdat.Mass,sizeof(float),Nmass,fp);
	SKIP2;
	CHECK;
		
	fclose(fp);
	//fprintf(logfile,"DM particles loaded\n");
	#undef SKIP
	#undef SKIP2
	#undef CHECK
}

void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int Nvir[3],i,j;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,14*4L,SEEK_CUR);
		fread(Nvir,sizeof(int),3,fp);
		for(j=0;j<3;j++)
		Mvir[i][j]=Nvir[j]*partmass;
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}
