#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"

#define NBIN_MAX 60
#define NBIN_MIN 5
#define Factor_RMIN 1e-2
#define Factor_relax 3

#define NDIV 200
LINKLIST ll;
#ifdef PERIODIC_BDR
#define MAKELL(pos,np,ndiv) make_linklist_box(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_shift_gridid(x,NDIV)
#else
#define MAKELL(pos,np,ndiv) make_linklist(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_round_gridid(x,NDIV)
#endif

CATALOGUE Cat;
GASHALOCAT GCat;

typedef struct 
{
	int nbin;
	float rmax;
	float rmin;//min radius of radial log-bins; equals to max(Factor_RMIN*rmax,SofteningHalo); note this is not saved.
	float CoM[3];//halo CoM,for comparison
	float Cen[3];//Center of mainsub,the profiles are w.r.t this center
	int flag_fakehalo;//set to 1 when halo is not self-bound,to 0 otherwise.
	
	float (*vm_xyz)[3];//(not physical)
	float (*vd_xyz)[3];
	float (*vm_rtf)[3];//spherical coordinate
	float (*vd_rtf)[3];
	float *Um;
	float (*gvm_xyz)[3];
	float (*gvd_xyz)[3];
	float (*gvm_rtf)[3];//(mean of vel)
	float (*gvd_rtf)[3];//velocity dispersion
	
} DYNAMPROFILE;  //profile for all particles in radial bin (this+back+other)

void dynam_halo_bin(int grpid,DYNAMPROFILE *haloprof);
void dynam_ghalo_bin(int grpid,DYNAMPROFILE *haloprof);
void free_dynamprof(DYNAMPROFILE *haloprof);
void load_halos_bin(DYNAMPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir);
void save_dynam_prof(DYNAMPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir);

int main(int argc,char **argv)
{
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	char gasdir[1024];
	
int Nsnap=99,i;
int grpid;
DYNAMPROFILE *haloprof;

Nsnap=atoi(argv[1]);

logfile=stdout;
sprintf(outputdir,"%s/profile/logbin",SUBCAT_DIR);	

	load_particle_data(Nsnap,SNAPSHOT_DIR);	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);

	#ifndef GRPINPUT_INDEX
	fresh_ID2Index(&Cat,-1); 	//fofcat of JING's data are originally PIndex rather than PID
	#endif
	prepare_ind2halo(&Cat);

	MAKELL(Pdat.Pos,NP_DM,NDIV);
	printf("finished linkedlist\n");

haloprof=mymalloc(sizeof(DYNAMPROFILE)*Cat.Ngroups);
load_halos_bin(haloprof,Cat.Ngroups,Nsnap,outputdir);
#pragma omp parallel for
for(grpid=0;grpid<Cat.Ngroups;grpid++)
{	
	if(haloprof[grpid].nbin)
	dynam_halo_bin(grpid,haloprof+grpid);
}
free_catalogue(&Cat);
	free_linklist(&ll);

sprintf(gasdir,"%s/gascat",SUBCAT_DIR);
load_gas_data(Nsnap,SNAPSHOT_DIR);
load_gashalocat(Nsnap,&GCat,GASCAT_DIR);
	#ifndef GRPINPUT_INDEX
	fresh_gasID2Index(&GCat,-1); 	//fofcat of JING's data are originally PIndex rather than PID
	#endif
	GCat.ID2Halo=mymalloc(sizeof(int)*NP_GAS);
	prepare_gasind2halo(&GCat);

	MAKELL(Gdat.Pos,NP_GAS,NDIV);
	printf("finished linkedlist\n");
#pragma omp parallel for
for(grpid=0;grpid<GCat.Ngroups;grpid++)
{	
	if(haloprof[grpid].nbin)
	dynam_ghalo_bin(grpid,haloprof+grpid);
}
	free_linklist(&ll);
free(GCat.ID2Halo);
free_gashalocat(&GCat);
save_dynam_prof(haloprof,GCat.Ngroups,Nsnap,outputdir);
for(grpid=0;grpid<Cat.Ngroups;grpid++)
	free_dynamprof(haloprof+grpid);
free(haloprof);

return 0;
}

void dynam_halo_bin(int grpid,DYNAMPROFILE *haloprof)
{
	float *cen,rmax,rmin,dr,vr,vt,vf;
	float (*vm_xyz)[3],(*vsm_xyz)[3],(*vm_rtf)[3],(*vsm_rtf)[3];
    double dx[3],er[3],et[3],ef[3];//unit vectors of spherical coord.
		  //~ *Um_this,*Um_other,*Um_back;
	int    *np,nbin,dr_bin;
	double span_order;

	int i,j,k,l,pid;
	int subbox_grid[3][2];

	if(haloprof->flag_fakehalo)
		cen=haloprof->CoM;
	else
		cen=haloprof->Cen;
		
	rmax=haloprof->rmax;
	rmin=haloprof->rmin;
	span_order=log(rmax/rmin);	
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rmax-ll.range[i][0])/ll.step[i]);
	subbox_grid[i][1]=floor((cen[i]+rmax-ll.range[i][0])/ll.step[i]);
	}
	nbin=haloprof->nbin;
	vm_xyz=mymalloc(sizeof(float)*nbin*3);
	vsm_xyz=mymalloc(sizeof(float)*nbin*3);
	vm_rtf=mymalloc(sizeof(float)*nbin*3);
	vsm_rtf=mymalloc(sizeof(float)*nbin*3);
	np=mymalloc(sizeof(float)*nbin);
	for(i=0;i<nbin;i++)
	{
		vm_xyz[i][0]=vm_xyz[i][1]=vm_xyz[i][2]=0.;
		vsm_xyz[i][0]=vsm_xyz[i][1]=vsm_xyz[i][2]=0.;
		np[i]=0;
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc(&ll,FIXGRID(i),FIXGRID(j),FIXGRID(k));
				while(pid>=0)
				{
					dr=distance(Pdat.Pos[pid],cen);
					if(dr<rmax&&dr>rmin)
					{
						dr_bin=floor(log(dr/rmin)/span_order*nbin);
						if(dr_bin>=nbin)	dr_bin=nbin-1;
						if(dr_bin<0)	dr_bin=0;
						for(l=0;l<3;l++)
							{
								vm_xyz[dr_bin][l]+=Pdat.Vel[pid][l];
								vsm_xyz[dr_bin][l]+=Pdat.Vel[pid][l]*Pdat.Vel[pid][l];
							}
							np[dr_bin]++;
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<nbin;i++)
		for(j=0;j<3;j++)
		{
		vm_xyz[i][j]/=np[i];
		vsm_xyz[i][j]/=np[i];
		vsm_xyz[i][j]-=vm_xyz[i][j]*vm_xyz[i][j];
		}		
	for(i=0;i<nbin;i++)
	{
		vm_rtf[i][0]=vm_rtf[i][1]=vm_rtf[i][2]=0.;
		vsm_rtf[i][0]=vsm_rtf[i][1]=vsm_rtf[i][2]=0.;
		np[i]=0;
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc(&ll,FIXGRID(i),FIXGRID(j),FIXGRID(k));
				while(pid>=0)
				{
					dr=distance(Pdat.Pos[pid],cen);
					if(dr<rmax&&dr>rmin)
					{
						dr_bin=floor(log(dr/rmin)/span_order*nbin);
						if(dr_bin>=nbin)	dr_bin=nbin-1;
						if(dr_bin<0)	dr_bin=0;
						for(l=0;l<3;l++) dx[l]=Pdat.Pos[pid][l]-cen[l];
						spherical_basisD(dx,er,et,ef);
						vr=vt=vf=0.;
						for(l=0;l<3;l++)
							{
								vr+=(Pdat.Vel[pid][l]-vm_xyz[dr_bin][l])*er[l];
								vt+=(Pdat.Vel[pid][l]-vm_xyz[dr_bin][l])*et[l];
								vf+=(Pdat.Vel[pid][l]-vm_xyz[dr_bin][l])*ef[l];
							}
								vm_rtf[dr_bin][0]+=vr;
								vm_rtf[dr_bin][1]+=vt;
								vm_rtf[dr_bin][2]+=vf;
								vsm_rtf[dr_bin][0]+=vr*vr;
								vsm_rtf[dr_bin][1]+=vt*vt;
								vsm_rtf[dr_bin][2]+=vf*vf;
							np[dr_bin]++;
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<nbin;i++)
		for(j=0;j<3;j++)
		{
		vm_rtf[i][j]/=np[i];
		vsm_rtf[i][j]/=np[i];
		vsm_rtf[i][j]-=vm_rtf[i][j]*vm_rtf[i][j];
		}
	haloprof->vm_xyz=vm_xyz;
	haloprof->vm_rtf=vm_rtf;
	haloprof->vd_xyz=vsm_xyz;
	haloprof->vd_rtf=vsm_rtf;
	free(np);
}

void dynam_ghalo_bin(int grpid,DYNAMPROFILE *haloprof)
{
	float *cen,rmax,rmin,dr,vr,vt,vf;
	float (*vm_xyz)[3],(*vm_rtf)[3],
		  (*vsm_xyz)[3],(*vsm_rtf)[3];
	double dx[3],er[3],et[3],ef[3];//unit vectors of spherical coord.
	float *Um;
	int    *np,nbin,dr_bin;
	double span_order;

	int i,j,k,l,pid;
	int subbox_grid[3][2];

	if(haloprof->flag_fakehalo)
		cen=haloprof->CoM;
	else
		cen=haloprof->Cen;
		
	rmax=haloprof->rmax;
	rmin=haloprof->rmin;
	span_order=log(rmax/rmin);	
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rmax-ll.range[i][0])/ll.step[i]);
	subbox_grid[i][1]=floor((cen[i]+rmax-ll.range[i][0])/ll.step[i]);
	}
	nbin=haloprof->nbin;
	vm_xyz=mymalloc(sizeof(float)*nbin*3);
	vsm_xyz=mymalloc(sizeof(float)*nbin*3);
	vm_rtf=mymalloc(sizeof(float)*nbin*3);
	vsm_rtf=mymalloc(sizeof(float)*nbin*3);
	Um=mymalloc(sizeof(float)*nbin);
	np=mymalloc(sizeof(float)*nbin);
	for(i=0;i<nbin;i++)
	{
		vm_xyz[i][0]=vm_xyz[i][1]=vm_xyz[i][2]=0.;
		vsm_xyz[i][0]=vsm_xyz[i][1]=vsm_xyz[i][2]=0.;
		Um[i]=0.;
		np[i]=0;
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc(&ll,FIXGRID(i),FIXGRID(j),FIXGRID(k));
				while(pid>=0)
				{
					dr=distance(Gdat.Pos[pid],cen);
					if(dr<rmax&&dr>rmin)
					{
						dr_bin=floor(log(dr/rmin)/span_order*nbin);
						if(dr_bin>=nbin)	dr_bin=nbin-1;
						if(dr_bin<0)	dr_bin=0;
						for(l=0;l<3;l++)
							{
								vm_xyz[dr_bin][l]+=Gdat.Vel[pid][l];
								vsm_xyz[dr_bin][l]+=Gdat.Vel[pid][l]*Gdat.Vel[pid][l];
							}
							Um[dr_bin]+=Gdat.U[pid];
							np[dr_bin]++;
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<nbin;i++)
	{
		for(j=0;j<3;j++)
		{
		vm_xyz[i][j]/=np[i];
		vsm_xyz[i][j]/=np[i];
		vsm_xyz[i][j]-=vm_xyz[i][j]*vm_xyz[i][j];
		}
		Um[i]/=np[i];
	}			
	for(i=0;i<nbin;i++)
	{
		vm_rtf[i][0]=vm_rtf[i][1]=vm_rtf[i][2]=0.;
		vsm_rtf[i][0]=vsm_rtf[i][1]=vsm_rtf[i][2]=0.;
		np[i]=0;
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc(&ll,FIXGRID(i),FIXGRID(j),FIXGRID(k));
				while(pid>=0)
				{
					dr=distance(Gdat.Pos[pid],cen);
					if(dr<rmax&&dr>rmin)
					{
						dr_bin=floor(log(dr/rmin)/span_order*nbin);
						if(dr_bin>=nbin)	dr_bin=nbin-1;
						if(dr_bin<0)	dr_bin=0;
						for(l=0;l<3;l++) dx[l]=Gdat.Pos[pid][l]-cen[l];
						spherical_basisD(dx,er,et,ef);
						vr=vt=vf=0.;
						for(l=0;l<3;l++)
							{
								vr+=(Gdat.Vel[pid][l]-vm_xyz[dr_bin][l])*er[l];
								vt+=(Gdat.Vel[pid][l]-vm_xyz[dr_bin][l])*et[l];
								vf+=(Gdat.Vel[pid][l]-vm_xyz[dr_bin][l])*ef[l];
							}
								vm_rtf[dr_bin][0]+=vr;
								vm_rtf[dr_bin][1]+=vt;
								vm_rtf[dr_bin][2]+=vf;
								vsm_rtf[dr_bin][0]+=vr*vr;
								vsm_rtf[dr_bin][1]+=vt*vt;
								vsm_rtf[dr_bin][2]+=vf*vf;
							np[dr_bin]++;
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<nbin;i++)
		for(j=0;j<3;j++)
		{
		vm_rtf[i][j]/=np[i];
		vsm_rtf[i][j]/=np[i];
		vsm_rtf[i][j]-=vm_rtf[i][j]*vm_rtf[i][j];
		}
	haloprof->gvm_xyz=vm_xyz;
	haloprof->gvm_rtf=vm_rtf;
	haloprof->gvd_xyz=vsm_xyz;
	haloprof->gvd_rtf=vsm_rtf;
	haloprof->Um=Um;
	free(np);
}

void free_dynamprof(DYNAMPROFILE *haloprof)
{
	if(haloprof->nbin)
	{
		free(haloprof->vm_xyz);
		free(haloprof->vd_xyz);
		free(haloprof->vm_rtf);
		free(haloprof->vd_rtf);
		free(haloprof->Um);
		free(haloprof->gvm_xyz);
		free(haloprof->gvd_xyz);
		free(haloprof->gvm_rtf);
		free(haloprof->gvd_rtf);
	}
}

void load_halos_bin(DYNAMPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/halo_size_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Nsubs;i++)
	{
		fread(&haloprof[i].nbin,sizeof(int),1,fp);
		fread(&haloprof[i].rmax,sizeof(float),1,fp);
		haloprof[i].rmin=haloprof[i].rmax*Factor_RMIN;
		if(haloprof[i].rmin<SofteningHalo)
			haloprof[i].rmin=SofteningHalo;
		fseek(fp,sizeof(float)*5,SEEK_CUR);
		fread(haloprof[i].CoM,sizeof(float),3,fp);
		fread(haloprof[i].Cen,sizeof(float),3,fp);
		fseek(fp,sizeof(float)*10,SEEK_CUR);
		fread(&haloprof[i].flag_fakehalo,sizeof(int),1,fp);
	}
	fclose(fp);
}

void save_dynam_prof(DYNAMPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i,nbin;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/dynam_prof_%03d.sph",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Nsubs;i++)
	{
		nbin=haloprof[i].nbin;
		fwrite(haloprof[i].vm_xyz,sizeof(float),nbin*3,fp);
		fwrite(haloprof[i].vd_xyz,sizeof(float),nbin*3,fp);
		fwrite(haloprof[i].vm_rtf,sizeof(float),nbin*3,fp);
		fwrite(haloprof[i].vd_rtf,sizeof(float),nbin*3,fp);
		fwrite(haloprof[i].Um,sizeof(float),nbin,fp);
		fwrite(haloprof[i].gvm_xyz,sizeof(float),nbin*3,fp);
		fwrite(haloprof[i].gvd_xyz,sizeof(float),nbin*3,fp);
		fwrite(haloprof[i].gvm_rtf,sizeof(float),nbin*3,fp);
		fwrite(haloprof[i].gvd_rtf,sizeof(float),nbin*3,fp);
	}
	fclose(fp);
}
