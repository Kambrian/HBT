//output data to check for the Mass-concentration relation
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);
int load_halo_concentration(float *halocon,int Nsnap)
{
	int grpid,Ngroups,Ngroups2,halostatus;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_param_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	fread(&Ngroups,sizeof(int),1,fp);
	for(grpid=0;grpid<Ngroups;grpid++)
	{
		fread(&halostatus,sizeof(int),1,fp);
		fseek(fp,4*sizeof(float),SEEK_CUR);
		fread(halocon+grpid,sizeof(float),1,fp);
		if(halostatus!=0) halocon[grpid]=-1; //set concentration to -1 if fitting not successful
		fseek(fp,sizeof(float)*2,SEEK_CUR);
	}	
	fread(&Ngroups2,sizeof(int),1,fp);
	if(Ngroups2!=Ngroups)
	{
		printf("error:Ngroups=%d,%d do not match when loading \n %s\n" 
			"probably file corruption or different file format\n",
			Ngroups,Ngroups2,buf);
		exit(-1);
	}
	fclose(fp);
	return Ngroups;
}
int read_Ngroups(int Nsnap)
{
FILE *fd;
char buf[1024];
int Ngroups;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(int), 1, fd);
  
  fclose(fd);
  return Ngroups;
}
int main(int argc,char **argv)
{
	FILE *fp;
	char buf[1024];
	int Nsnap,i,Nmin,Nwrt,Ngroups;
	float *halocon,(*Mvir)[3],(*Rvir)[3];
	int Mmin;
	float MassMin,partmass;
	
	logfile=stdout;
	if(argc!=3)
	{
	printf("usage:%s [Snap] [HaloLenMin]\n",argv[0]);
	exit(1);
	}
	Nsnap=atoi(argv[1]);
	Mmin=atoi(argv[2]);
	
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	MassMin=Mmin*partmass;
	
	Ngroups=read_Ngroups(Nsnap);
	Mvir=mymalloc(sizeof(float)*3*Ngroups);
	Rvir=mymalloc(sizeof(float)*3*Ngroups);
	halocon=mymalloc(sizeof(float)*Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,Ngroups,Nsnap);
	load_halo_concentration(halocon,Nsnap);
	
	Nmin=0;
	for(i=0;i<Ngroups;i++)
	{
		if(Mvir[i][0]>MassMin&&halocon[i]>0)
		Nmin++;
	}
	printf("%d groups will be written\n",Nmin);
	sprintf(buf,"%s/anal/M_c_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(&Nmin,sizeof(int),1,fp);
	Nwrt=0;
	for(i=0;i<Ngroups;i++)
	{
		if(Mvir[i][0]>MassMin&&halocon[i]>0)
		{
		fwrite(&Mvir[i][0],sizeof(float),1,fp);
		fwrite(halocon+i,sizeof(float),1,fp);
		Nwrt++;
		}
	}
	fwrite(&Nwrt,sizeof(int),1,fp);
	if(Nwrt!=Nmin)
	{
		printf("error: Ngroups written=%d, %d\n",Nwrt,Nmin);
		exit(1);
	}
	fclose(fp);
	myfree(Mvir);
	myfree(Rvir);
	myfree(halocon);
	
	return 0;
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
