#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define _TIDAL  //calculate using Rtidal as scale, otherwise use virial radius
#define _COM //calculate using (core-)CoM as sub position. otherwise use most bound
#define NDIV 256
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
int Nsnap=99,i,j,k,l,pid,subid=1,grpid=0;

int Ngas,Ndm,Nother,Nmass,dummy,dummy2;
float *cen,rvir,r_in,r_out,r_mid,L_in,L_mid,L_out,V_in,V_mid,V_out,f_in,f_out,f_mid;
//~ int mapxy[NGRID][NGRID]={0},mapxz[NGRID][NGRID]={0},mapyz[NGRID][NGRID]={0};
int *PIndex,grid[3];
//~ double com[3];

int subbox_grid[3][2];
float range[3][2],step[3],dr;

float lowerlim,upperlim,Lx;
double rfun(double,void *);
double integ,err;
gsl_function F;
size_t neval;

float *Rtidal;
char buf[1024];
FILE * fpsnap,*fpsublx,*fpsubraw,*fpsubvel;
logfile=stdout;

F.function=&rfun;
	
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	#ifdef _TIDAL
	Rtidal=mymalloc(sizeof(float)*SubCat.Nsubs);
	load_tidal_radius(Nsnap,Rtidal,SubCat.Nsubs,inputdir);
	#endif
	load_particle_data(Nsnap,snapdir);fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);
	PIndex=Cat.PIDorIndex+Cat.Offset[grpid];
	//~ for(i=0;i<3;i++)
		//~ for(j=0;j<2;j++)
			//~ range[i][j]=Pdat.Pos[PIndex[0]][i];
	//~ for(i=1;i<Cat.Len[grpid];i++)
	//~ {
		//~ pid=PIndex[i];
		//~ for(j=0;j<3;j++)
		//~ {
		//~ if(Pdat.Pos[pid][j]<range[j][0])
			//~ range[j][0]=Pdat.Pos[pid][j];
		//~ else if(Pdat.Pos[pid][j]>range[j][1])
			//~ range[j][1]=Pdat.Pos[pid][j];
		//~ }
	//~ }


	//~ printf("size: %f,%f\n%f,%f\n%f,%f\n",range[0][1],range[0][0],range[1][1],range[1][0],range[2][1],range[2][0]);
	
	sprintf(buf,"%s/snapshot_%03d",snapdir,Nsnap);myfopen(fpsnap,buf,"r");
	sprintf(buf,"%s/subLx_%03d",outputdir,Nsnap);myfopen(fpsublx,buf,"w");
	sprintf(buf,"%s/subLxraw_%03d",outputdir,Nsnap);myfopen(fpsubraw,buf,"w");
	sprintf(buf,"%s/subvel_%03d",outputdir,Nsnap);myfopen(fpsubvel,buf,"w");

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
	f_in=1;f_out=1.5;
	f_mid=pow(f_in*f_in*f_in+(f_out*f_out*f_out-f_in*f_in*f_in)/2,1.0/3);
	
	MAKELL(gPos,Ngas,NDIV);
	//printf("%d,%d,%g\n",hoc[10][9][133],hoc[77][2][1],Hubble);
	printf("finished linkedlist\n");
	
for(grpid=0;grpid<SubCat.Ngroups;grpid++)
{
for(subid=SubCat.GrpOffset_Sub[grpid]+1;subid<SubCat.GrpOffset_Sub[grpid+1];subid++)
{
	#ifdef _COM
	cen=SubCat.CoM[subid];
	#else
	cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
	#endif
	#ifdef _TIDAL
	rvir=Rtidal[subid];
	#else
	rvir=pow(G*SubCat.SubLen[subid]*headerA.mass[1]/100/headerA_Hz/headerA_Hz,1.0/3);
	#endif
	r_in=f_in*rvir;	r_out=f_out*rvir;	r_mid=f_mid*rvir;
	L_in=0;V_in=0;	L_mid=0;V_mid=0;	L_out=0;V_out=0;
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-r_out-ll.range[i][0])/ll.step[i]);
	subbox_grid[i][1]=floor((cen[i]+r_out-ll.range[i][0])/ll.step[i]);
	}
		for(l=subbox_grid[0][0];l<subbox_grid[0][1]+1;l++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				i=linklist_get_hoc(&ll,FIXGRID(l),FIXGRID(j),FIXGRID(k));
				while(i>=0)
				{
					dr=distance(gPos[i],cen);
						if(dr<r_out)
						{
							if(gU[i]<5e3)//[lower,upper]~[20,80],the integral is 10 order of magnitude lower than that in [1,4],so omit these non-X-ray gas (T<2e5K or 0.02keV)
								Lx=0;
							else
							{
								lowerlim=0.5e6/3.48/gU[i]; //0.5~2kev
								upperlim=2e6/3.48/gU[i];
								gsl_integration_qng(&F,lowerlim,upperlim,0,0.001,&integ,&err,&neval);
								Lx=headerA.mass[0]*grho[i]/(headerA.time*headerA.time*headerA.time)*sqrt(gU[i])*integ*3040.7;
							}
							if(dr<r_in)//innermost
							{
									L_in+=Lx;
									V_in+=1./grho[i];
							}
							else if(dr<r_mid)//middle belt
							{
									L_mid+=Lx;
									V_mid+=1./grho[i];
							}
								else//outer belt
							{
									L_out+=Lx;
									V_out+=1./grho[i];
							}
						}
						i=ll.list[i];
				}
			}
fprintf(fpsublx,"%d,%g,%g,%g\n",subid,L_in-L_mid*V_in/V_mid,L_in-L_out*V_in/V_out,L_in-(L_mid+L_out)*V_in/(V_mid+V_out));
fprintf(fpsubraw,"%d,%g,%g,%g,%g,%g,%g\n",subid,L_in,L_mid,L_out,V_in,V_mid,V_out);
fprintf(fpsubvel,"%d,%g,%g,%g\n",subid,Pdat.Vel[SubCat.PSubArr[subid][0]][0],Pdat.Vel[SubCat.PSubArr[subid][0]][1],Pdat.Vel[SubCat.PSubArr[subid][0]][2]);
}
}

fclose(fpsublx);
fclose(fpsubraw);
fclose(fpsubvel);

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
