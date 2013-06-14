#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define _TIDAL   //define this to use Rtidal as scale,otherwise use virial radius
#define _COM     //define this to use (core-)CoM as sub center, otherwise use most bound
#define NDIV 200
#define NBIN 10
LINKLIST ll;
#ifdef PERIODIC_BDR
#define MAKELL(pos,np,ndiv) make_linklist_box(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_shift_gridid(x,NDIV)
#else
#define MAKELL(pos,np,ndiv) make_linklist(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_round_gridid(x,NDIV)
#endif
void load_tidal_radius(int Nsnap,float *rtidal, int Nsubs,char*tidaldir);

int main()
{
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	#ifdef _TIDAL
			sprintf(outputdir,"%s/anal/LxTidal",SUBCAT_DIR);	
		#else
			sprintf(outputdir,"%s/anal/LxVir",SUBCAT_DIR);	
	#endif
	#ifdef _COM
			sprintf(outputdir,"%s/CoM",outputdir);
	#endif

CATALOGUE Cat;
SUBCATALOGUE SubCat;
int *gPID;
float (*gPos)[3],(*gVel)[3],*gU,*grho;
int Nsnap=99,i,j,k,pid,subid,grpid;

int Ngas,Ndm,Nother,Nmass,dummy,dummy2;
float *cen,rvir,L_div[NBIN],Rgho_div[NBIN],T_div[NBIN],r_div[NBIN];
int n_div[NBIN];

int subbox_grid[3][2],dr_bin;
float range[3][2],step[3],dr;
float *Rtidal;

float lowerlim,upperlim,Lx;
double rfun(double,void *);
double integ,err;
gsl_function F;
size_t neval;

char buf[1024];
FILE * fpsnap,*fpsubprof;
logfile=stdout;

F.function=&rfun;
	
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	#ifdef _TIDAL
	Rtidal=mymalloc(sizeof(float)*SubCat.Nsubs);
	load_tidal_radius(Nsnap,Rtidal,SubCat.Nsubs,inputdir);
	#endif
	load_particle_data(Nsnap,snapdir);fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);

	sprintf(buf,"%s/snapshot_%03d",snapdir,Nsnap);myfopen(fpsnap,buf,"r");
	sprintf(buf,"%s/subprof_%03d",outputdir,Nsnap);myfopen(fpsubprof,buf,"w");

#define _SKIP fread(&dummy,sizeof(dummy),1,fpsnap)
#define _SKIP2 fread(&dummy2,sizeof(dummy2),1,fpsnap)
#define _CHECK if(dummy!=dummy2){printf("error!record brackets not match for Snap%d!\t%d,%d\n",Nsnap,dummy,dummy2);exit(1);} 	

	_SKIP;	fread(&headerA,sizeof(headerA),1,fpsnap);	_SKIP2;	_CHECK;
	Ngas=headerA.npart[0];	Ndm=headerA.npart[1];	Nother=headerA.npart[2]+headerA.npart[3]+headerA.npart[4]+headerA.npart[5];
	for(Nmass=0,i=0;i<6;i++)
	{
		if(!headerA.mass[i]) Nmass+=headerA.npart[i];
	}
	printf("%f,%f,%f,%f,%f,%f,%d,%d\n",headerA.mass[0],headerA.mass[1],headerA.mass[2],headerA.mass[3],headerA.mass[4],headerA.mass[5],Nmass,Nother);
	gPID=malloc(sizeof(int)*Ngas);
	gPos=malloc(sizeof(float)*3*Ngas);
	gVel=malloc(sizeof(float)*3*Ngas);
	gU=malloc(sizeof(float)*Ngas);
	grho=malloc(sizeof(float)*Ngas);
	
	_SKIP;
	fread(gPos,sizeof(float)*3,Ngas,fpsnap);
	fseek(fpsnap,(Ndm+Nother)*sizeof(float)*3,SEEK_CUR);
	_SKIP2;
	_CHECK;
	
	_SKIP;
	fread(gVel,sizeof(float)*3,Ngas,fpsnap);
	fseek(fpsnap,(Ndm+Nother)*sizeof(float)*3,SEEK_CUR);
	_SKIP2;
	_CHECK;

	_SKIP;
	fread(gPID,sizeof(int),Ngas,fpsnap);
	fseek(fpsnap,(Ndm+Nother)*sizeof(int),SEEK_CUR);
	_SKIP2;
	_CHECK;
	
	if(Nmass)
	{
		_SKIP;
		fseek(fpsnap,Nmass*sizeof(float),SEEK_CUR);
		_SKIP2;
		_CHECK;
	}
	
	_SKIP;
	fread(gU,sizeof(float),Ngas,fpsnap);
	_SKIP2;
	_CHECK;

	_SKIP;
	fread(grho,sizeof(float),Ngas,fpsnap);
	_SKIP2;
	_CHECK;
	fclose(fpsnap);

	//~ for(i=0;i<10;i++)	
		//~ f_div[i]=0.1*(i+1);

	MAKELL(gPos,Ngas,NDIV);
	//printf("%d,%d,%g\n",hoc[100][100][110],hoc[77][92][101],Hubble);
	printf("finished linkedlist\n");
	//~ printf("%g,%g,%g\n",gPos[0][0],gPos[0][1],gPos[0][2]);
for(grpid=0;grpid<SubCat.Ngroups;grpid++)
{
for(subid=SubCat.GrpOffset_Sub[grpid]+1;subid<SubCat.GrpOffset_Sub[grpid+1];subid++)
{
	for(i=0;i<NBIN;i++)	
	{
		L_div[i]=0;
		T_div[i]=0;
		Rgho_div[i]=0;
		n_div[i]=0;
		r_div[i]=0;
	}
	#ifdef _COM
	cen=SubCat.CoM[subid];
	#else
	cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
	#endif
	#ifdef _TIDAL
	rvir=Rtidal[subid];
	#else
	rvir=pow(G*SubCat.SubLen[subid]*header.mass[1]/100/header.Hz/header.Hz,1.0/3)/header.time;
	#endif
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rvir-ll.range[i][0])/ll.step[i]);
	subbox_grid[i][1]=floor((cen[i]+rvir-ll.range[i][0])/ll.step[i]);
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc(&ll,FIXGRID(i),FIXGRID(j),FIXGRID(k));
				while(pid>=0)
				{
					dr=distance(gPos[pid],cen);
					if(dr<rvir)
					{
						dr_bin=floor(dr/rvir*NBIN);
						if(gU[pid]<5e3)//[lower,upper]~[20,80],the integral is 10 order of magnitude lower than that in [1,4],so omit these non-X-ray gas (T<2e5K)
							Lx=0;
						else
						{
							lowerlim=0.5e6/3.48/gU[pid]; //0.5~2kev
							upperlim=2e6/3.48/gU[pid];
							gsl_integration_qng(&F,lowerlim,upperlim,0,0.001,&integ,&err,&neval);
							Lx=headerA.mass[0]*grho[pid]/(headerA.time*headerA.time*headerA.time)*sqrt(gU[pid])*integ*3040.7;
						}
						L_div[dr_bin]+=Lx;
						Rgho_div[dr_bin]+=grho[pid];
						T_div[dr_bin]+=gU[pid];
						r_div[dr_bin]+=dr;
						n_div[dr_bin]++;
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<NBIN;i++)
	{
		if(n_div[i])
		{
		L_div[i]/=n_div[i];
		Rgho_div[i]/=n_div[i];
		T_div[i]=40.38*T_div[i]/n_div[i];
		r_div[i]=r_div[i]/n_div[i]/rvir;
		}
		else
		{
			r_div[i]=(i+0.5)/NBIN;
		}
	}
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",L_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",Rgho_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",T_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",r_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%d\t",n_div[i]);
		fprintf(fpsubprof,"%g\t%d\t%d\n",rvir,SubCat.SubLen[subid],subid);
		//~ fflush(fpsubprof);
		//~ printf("%d\n",subid);
}
}
		
fclose(fpsubprof);
return 0;
}

double rfun(double x,void * param)
{
	return 1.0/(pow(x,0.4)*exp(x));
}


//~ float rfun(float x)
//~ {
	//~ return 1.0/(pow(x,0.4)*exp(x));
//~ }
void load_tidal_radius(int Nsnap,float *rtidal, int Nsubs,char*tidaldir)
{
	char buf[1024];
FILE *fp;
	int i;
	
	sprintf(buf,"%s/tidal_radius_%03d",tidaldir,Nsnap);
	if(!(fp=fopen(buf,"r")))
	{
		fprintf(logfile,"Error opening file '%s'\n",buf);
		exit(1);
	}
	for(i=0;i<Nsubs;i++)	
	{
		fscanf(fp,"%e\n",rtidal+i);
		if(ferror(fp))
		{
		fprintf(logfile,"error reading file %s: i=%d,Nsubs=%d\n",buf,i,Nsubs);
		exit(1);
		}
	}
	fclose(fp);
}
