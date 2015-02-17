//halo abundance (number of subhalos inside virial radius above some (relative) mass-threshold)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

//~ #define NORM  //produce normalized massfunction (in terms Msub/Mhost rather than Msub)
#define MINSUBMASS 1e-4
#define RMIN 0
#define RMAX 1	 //statistics done in RMIN*rvi<r<RMAX*rvir


double partmass;
float (*Rvir)[3];
int (*Nvir)[3];
CATALOGUE Cat;
SUBCATALOGUE SubCat;
int VirType;

void logspace(double xmin,double xmax,int N,float *x);
void load_halo_virial_size(int Nvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);
void makell_sub();
int collect_subnumber(int grpid,float rmin,float rmax);
void freell_sub();
void load_halos_param(float (**p2halopars)[7],int **p2halostatus,int Nsnap);
int main(int argc,char **argv)
{
	int Nsnap,i,j;
	float MvirMin;
	float (*halopars)[7];
	int *halostatus;

	FILE *fp;
	char buf[1024];
	
	logfile=stdout;
	if(argc!=2)
	{
	printf("usage:%s [Snap]\n",argv[0]);
	exit(1);
	}
	Nsnap=atoi(argv[1]);
	
	VirType=0;
		
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/massfun",SUBCAT_DIR);
	mkdir(outputdir,0755);
	
	#ifdef NORM
	sprintf(buf,"%s/anal/massfun/abundanceN_%03d.%e",SUBCAT_DIR,Nsnap,MINSUBMASS);
	#else
	sprintf(buf,"%s/anal/massfun/abundance_%03d",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"w");

	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	Nvir=mymalloc(sizeof(int)*3*Cat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*Cat.Ngroups);
	load_halo_virial_size(Nvir,Rvir,(float)partmass,Cat.Ngroups,Nsnap);
	load_halos_param(&halopars,&halostatus,Nsnap);
	
	MvirMin=NBOUNDMIN/MINSUBMASS;

	makell_sub();
	int Ns,grpid;
	float Mvir;
	for(grpid=0;grpid<SubCat.Ngroups;grpid++)
	{
		#ifdef NORM
		if(Nvir[grpid][VirType]>MvirMin)
		#endif
		{
			Ns=collect_subnumber(grpid,RMIN*Rvir[grpid][VirType],RMAX*Rvir[grpid][VirType]);
			fwrite(&grpid,sizeof(int),1,fp);
			fwrite(&Ns,sizeof(int),1,fp);
			fwrite(&halostatus[grpid],sizeof(int),1,fp);
			fwrite(&halopars[grpid][0],sizeof(float),1,fp);
			fwrite(&halopars[grpid][2],sizeof(float),1,fp);
			fwrite(&halopars[grpid][4],sizeof(float),1,fp);
			fwrite(&halopars[grpid][6],sizeof(float),1,fp);
			Mvir=Nvir[grpid][0]*partmass;
			fwrite(&Mvir,sizeof(float),1,fp);
			fwrite(&Rvir[grpid][0],sizeof(float),1,fp);	
		}
	}
	
	fclose(fp);
	
	freell_sub();
	myfree(Nvir);
	myfree(Rvir);
	
	return 0;
}

void box_subcom()
{//move subcom to be inside one box
int i,j;
for(i=0;i<SubCat.Nsubs;i++)
{
	for(j=0;j<3;j++)
		if(SubCat.Property[i].CoM[j]<0||SubCat.Property[i].CoM[j]>BOXSIZE)
			SubCat.Property[i].CoM[j]=SubCat.Property[i].CoM[j]-floor(SubCat.Property[i].CoM[j]/BOXSIZE)*BOXSIZE;
}	
}
#define NDIV 200
static int hoc[NDIV][NDIV][NDIV],*ll;
static float range[3][2], step[3];
void makell_sub()
{
	int i,j,grid[3],np;
	
	#define POS(i,j) SubCat.Property[i].CoM[j]
	np=SubCat.Nsubs;
	printf("creating linked list..\n");
	ll=mymalloc(sizeof(int)*np);
	/*determining enclosing cube*/
	#ifdef PERIODIC_BDR
	box_subcom();
	#endif
	for(i=0;i<3;i++) //make this to be inside one periodic box
	{
	range[i][0]=0;
	range[i][1]=BOXSIZE;
	}
	//~ for(i=0;i<3;i++)
		//~ for(j=0;j<2;j++)
			//~ range[i][j]=POS(0,i);
	//~ for(i=1;i<np;i++)
		//~ for(j=0;j<3;j++)
		//~ {
			//~ if(POS(i,j)<range[j][0])
				//~ range[j][0]=POS(i,j);
			//~ else if(POS(i,j)>range[j][1])
				//~ range[j][1]=POS(i,j);
		//~ }
	for(j=0;j<3;j++)
		step[j]=(range[j][1]-range[j][0])/NDIV;
	
	/*initialize hoc*/
	int *phoc=&(hoc[0][0][0]);
	for(i=0;i<NDIV*NDIV*NDIV;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((POS(i,j)-range[j][0])/step[j]);
			if(grid[j]<0) 
				grid[j]=0;
			else if(grid[j]>=NDIV)
				grid[j]=NDIV-1;
		}
		ll[i]=hoc[grid[0]][grid[1]][grid[2]];
		hoc[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
	#undef POS(i,j)
}
void freell_sub()
{
	myfree(ll);
}
#ifdef PERIODIC_BDR
#define BOXGRID(i) ((i)<0?(i)%NDIV+NDIV:((i)>=NDIV?(i)%NDIV:(i)))
#else
#define BOXGRID(i) (i)
#endif
int collect_subnumber(int grpid,float rmin,float rmax)
{
	int i,j,k,ii,jj,kk,cenid,pid,subbox_grid[3][2],len;
	float *cen,rscale,dr;
	int MsubMin;
	
	MsubMin=ceilf(MINSUBMASS*Nvir[grpid][VirType]);

	len=0;
	cenid=SubCat.GrpOffset_Sub[grpid];
	if(SubCat.SubLen[cenid])
	{
		cen=SubCat.Property[cenid].CoM;
		rscale=rmax;
		for(i=0;i<3;i++)
		{
		subbox_grid[i][0]=floor((cen[i]-rscale-range[i][0])/step[i]);
		subbox_grid[i][1]=floor((cen[i]+rscale-range[i][0])/step[i]);
		#ifndef PERIODIC_BDR
		if(subbox_grid[i][0]<0)subbox_grid[i][0]=0;
		if(subbox_grid[i][1]>=NDIV)subbox_grid[i][1]=NDIV-1;
		#endif
		}
		//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
		for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		{
			ii=BOXGRID(i);  //periodic boundaries
			for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			{
				jj=BOXGRID(j);
				for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
				{	
					kk=BOXGRID(k);
					pid=hoc[ii][jj][kk];
					while(pid>=0)
					{
						dr=distance(SubCat.Property[pid].CoM,cen);
						if(pid!=cenid&&SubCat.SubLen[pid]&&dr<rmax&&dr>rmin)
						{
							#ifdef NORM
							if(SubCat.SubLen[pid]>=MsubMin)
							#endif
							len++;
						}
						pid=ll[pid];
					}
				}
			}
		}
	}
	return len;
}

void load_halo_virial_size(int Nvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int i,j;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,14*4L,SEEK_CUR);
		fread(Nvir+i,sizeof(int),3,fp);
		//~ for(j=0;j<3;j++)
		//~ Mvir[i][j]=Nvir[j]*partmass;
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}

void load_halos_param(float (**p2halopars)[7],int **p2halostatus,int Nsnap)
{
	int grpid,Ngroups,Ngroups2;
	char buf[1024];
	float (*halopars)[7];
	int *halostatus;
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_param_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	fread(&Ngroups,sizeof(int),1,fp);
	halostatus=mymalloc(sizeof(int)*Ngroups);
	halopars=mymalloc(sizeof(float)*Ngroups*7);
	for(grpid=0;grpid<Ngroups;grpid++)
	{
		fread(halostatus+grpid,sizeof(int),1,fp);
		fread(halopars[grpid],sizeof(float),7,fp);
	}	
	fread(&Ngroups2,sizeof(int),1,fp);
	fclose(fp);
	if(Ngroups!=Ngroups)
	{
		printf("error loading halo params: check failed %d, %d\n",Ngroups,Ngroups2);
		exit(1);
	}
	*p2halostatus=halostatus;
	*p2halopars=halopars;
}

void logspace(double xmin,double xmax,int N,float *x)
{
	int i;
	double dx;
	x[0]=xmin;x[N-1]=xmax;
	xmin=log(xmin);
	xmax=log(xmax);
	dx=exp((xmax-xmin)/(N-1));
	for(i=1;i<N-1;i++)
	{
		x[i]=x[i-1]*dx;
	}
}
