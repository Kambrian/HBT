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
	
	float *vm2_this;//(mean of vel)^2, need conversion to physical vel
	float *vm2_other;
	float *vm2_back;
	float *vrm2_this;//(mean of radial vel)^2
	float *vrm2_other;
	float *vrm2_back;
	//~ float (*vm_rtf)[3];//spherical coordinate
	float *vsm_this;//mean of square velocity
	float *vsm_other;
	float *vsm_back;
	float *vrsm_this;//mean of square radial velocity
	float *vrsm_other;
	float *vrsm_back;
	//~ float (*vd_rtf)[3];
	float *Um_this;
	float *Um_other;
	float *Um_back;
	float *gvm2_this;//(mean of vel)^2
	float *gvm2_other;
	float *gvm2_back;
	float *gvrm2_this;//(mean of radial vel)^2
	float *gvrm2_other;
	float *gvrm2_back;
	float *gvsm_this;//mean of square velocity
	float *gvsm_other;
	float *gvsm_back;
	float *gvrsm_this;//mean of square radial velocity
	float *gvrsm_other;
	float *gvrsm_back;
} DYNAMPROFILE;

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
	float *cen,rmax,rmin,dr;
	float (*vm_this)[3],(*vm_other)[3],(*vm_back)[3],
			*vm2_this,*vm2_other,*vm2_back,
		  *vsm_this,*vsm_other,*vsm_back;
		  //~ *Um_this,*Um_other,*Um_back;
	int    *np_this,*np_other,*np_back,nbin,dr_bin;
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
	vm2_this=mymalloc(sizeof(float)*nbin);
	vm2_other=mymalloc(sizeof(float)*nbin);
	vm2_back=mymalloc(sizeof(float)*nbin);
	vm_this=mymalloc(sizeof(float)*3*nbin);
	vm_other=mymalloc(sizeof(float)*3*nbin);
	vm_back=mymalloc(sizeof(float)*3*nbin);
	vsm_this=mymalloc(sizeof(float)*nbin);
	vsm_other=mymalloc(sizeof(float)*nbin);
	vsm_back=mymalloc(sizeof(float)*nbin);
	np_this=mymalloc(sizeof(float)*nbin);
	np_other=mymalloc(sizeof(float)*nbin);
	np_back=mymalloc(sizeof(float)*nbin);
	for(i=0;i<nbin;i++)
	{
		vm2_this[i]=vm2_other[i]=vm2_back[i]=0.;
		vm_this[i][0]=vm_other[i][0]=vm_back[i][0]=0.;
		vm_this[i][1]=vm_other[i][1]=vm_back[i][1]=0.;
		vm_this[i][2]=vm_other[i][2]=vm_back[i][2]=0.;
		vsm_this[i]=vsm_other[i]=vsm_back[i]=0.;
		np_this[i]=np_other[i]=np_back[i]=0;
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
						if(Cat.ID2Halo[pid]==grpid)
						{
							for(l=0;l<3;l++)
							{
								vm_this[dr_bin][l]+=Pdat.Vel[pid][l];
								vsm_this[dr_bin]+=Pdat.Vel[pid][l]*Pdat.Vel[pid][l];
							}
							np_this[dr_bin]++;
						}
						else if(Cat.ID2Halo[pid]==-1)
						{
							for(l=0;l<3;l++)
							{
								vm_back[dr_bin][l]+=Pdat.Vel[pid][l];
								vsm_back[dr_bin]+=Pdat.Vel[pid][l]*Pdat.Vel[pid][l];
							}
							np_back[dr_bin]++;
						}
						else
						{
							for(l=0;l<3;l++)
							{
								vm_other[dr_bin][l]+=Pdat.Vel[pid][l];
								vsm_other[dr_bin]+=Pdat.Vel[pid][l]*Pdat.Vel[pid][l];
							}
							np_other[dr_bin]++;
						}						
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<nbin;i++)
	{
		for(j=0;j<3;j++)
		{
		vm_this[i][j]/=np_this[i];
		vm_other[i][j]/=np_other[i];
		vm_back[i][j]/=np_back[i];
		vm2_this[i]+=vm_this[i][j]*vm_this[i][j];
		vm2_other[i]+=vm_other[i][j]*vm_other[i][j];
		vm2_back[i]+=vm_back[i][j]*vm_back[i][j];
		}
		vsm_this[i]/=np_this[i];
		vsm_other[i]/=np_other[i];
		vsm_back[i]/=np_back[i];
	}
	haloprof->vm2_this=vm2_this;
	haloprof->vm2_other=vm2_other;
	haloprof->vm2_back=vm2_back;
	haloprof->vsm_this=vsm_this;
	haloprof->vsm_other=vsm_other;
	haloprof->vsm_back=vsm_back;		
	free(vm_this);
	free(vm_other);
	free(vm_back);
	free(np_this);
	free(np_other);
	free(np_back);

}

void dynam_ghalo_bin(int grpid,DYNAMPROFILE *haloprof)
{
	float *cen,rmax,rmin,dr;
	float (*vm_this)[3],(*vm_other)[3],(*vm_back)[3],
			*vm2_this,*vm2_other,*vm2_back,
		  *vsm_this,*vsm_other,*vsm_back;
	float *Um_this,*Um_other,*Um_back;
	int    *np_this,*np_other,*np_back,nbin,dr_bin;
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
	vm2_this=mymalloc(sizeof(float)*nbin);
	vm2_other=mymalloc(sizeof(float)*nbin);
	vm2_back=mymalloc(sizeof(float)*nbin);
	vm_this=mymalloc(sizeof(float)*3*nbin);
	vm_other=mymalloc(sizeof(float)*3*nbin);
	vm_back=mymalloc(sizeof(float)*3*nbin);
	vsm_this=mymalloc(sizeof(float)*nbin);
	vsm_other=mymalloc(sizeof(float)*nbin);
	vsm_back=mymalloc(sizeof(float)*nbin);
	Um_this=mymalloc(sizeof(float)*nbin);
	Um_other=mymalloc(sizeof(float)*nbin);
	Um_back=mymalloc(sizeof(float)*nbin);
	np_this=mymalloc(sizeof(float)*nbin);
	np_other=mymalloc(sizeof(float)*nbin);
	np_back=mymalloc(sizeof(float)*nbin);
	for(i=0;i<nbin;i++)
	{
		vm2_this[i]=vm2_other[i]=vm2_back[i]=0.;
		vm_this[i][0]=vm_other[i][0]=vm_back[i][0]=0.;
		vm_this[i][1]=vm_other[i][1]=vm_back[i][1]=0.;
		vm_this[i][2]=vm_other[i][2]=vm_back[i][2]=0.;
		vsm_this[i]=vsm_other[i]=vsm_back[i]=0.;
		Um_this[i]=Um_other[i]=Um_back[i]=0.;
		np_this[i]=np_other[i]=np_back[i]=0;
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
						if(GCat.ID2Halo[pid]==grpid)
						{
							for(l=0;l<3;l++)
							{
								vm_this[dr_bin][l]+=Gdat.Vel[pid][l];
								vsm_this[dr_bin]+=Gdat.Vel[pid][l]*Gdat.Vel[pid][l];
							}
							Um_this[dr_bin]+=Gdat.U[pid];
							np_this[dr_bin]++;
						}
						else if(GCat.ID2Halo[pid]==-1)
						{
							for(l=0;l<3;l++)
							{
								vm_back[dr_bin][l]+=Gdat.Vel[pid][l];
								vsm_back[dr_bin]+=Gdat.Vel[pid][l]*Gdat.Vel[pid][l];
							}
							Um_back[dr_bin]+=Gdat.U[pid];
							np_back[dr_bin]++;
						}
						else
						{
							for(l=0;l<3;l++)
							{
								vm_other[dr_bin][l]+=Gdat.Vel[pid][l];
								vsm_other[dr_bin]+=Gdat.Vel[pid][l]*Gdat.Vel[pid][l];
							}
							Um_other[dr_bin]+=Gdat.U[pid];
							np_other[dr_bin]++;
						}						
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<nbin;i++)
	{
		for(j=0;j<3;j++)
		{
		vm_this[i][j]/=np_this[i];
		vm_other[i][j]/=np_other[i];
		vm_back[i][j]/=np_back[i];
		vm2_this[i]+=vm_this[i][j]*vm_this[i][j];
		vm2_other[i]+=vm_other[i][j]*vm_other[i][j];
		vm2_back[i]+=vm_back[i][j]*vm_back[i][j];
		}
		vsm_this[i]/=np_this[i];
		vsm_other[i]/=np_other[i];
		vsm_back[i]/=np_back[i];
		Um_this[i]/=np_this[i];
		Um_other[i]/=np_other[i];
		Um_back[i]/=np_back[i];
	}
	haloprof->Um_this=Um_this;
	haloprof->Um_other=Um_other;
	haloprof->Um_back=Um_back;
	haloprof->gvm2_this=vm2_this;
	haloprof->gvm2_other=vm2_other;
	haloprof->gvm2_back=vm2_back;
	haloprof->gvsm_this=vsm_this;
	haloprof->gvsm_other=vsm_other;
	haloprof->gvsm_back=vsm_back;		
	free(vm_this);
	free(vm_other);
	free(vm_back);
	free(np_this);
	free(np_other);
	free(np_back);		
}

void free_dynamprof(DYNAMPROFILE *haloprof)
{
	if(haloprof->nbin)
	{
		free(haloprof->vm2_this);
		free(haloprof->vm2_other);
		free(haloprof->vm2_back);
		free(haloprof->vsm_this);
		free(haloprof->vsm_other);
		free(haloprof->vsm_back);
		free(haloprof->Um_this);
		free(haloprof->Um_other);
		free(haloprof->Um_back);
		free(haloprof->gvm2_this);
		free(haloprof->gvm2_other);
		free(haloprof->gvm2_back);
		free(haloprof->gvsm_this);
		free(haloprof->gvsm_other);
		free(haloprof->gvsm_back);
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
	sprintf(buf,"%s/dynam_prof_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Nsubs;i++)
	{
		nbin=haloprof[i].nbin;
		fwrite(haloprof[i].vm2_this,sizeof(float),nbin,fp);
		fwrite(haloprof[i].vm2_other,sizeof(float),nbin,fp);
		fwrite(haloprof[i].vm2_back,sizeof(float),nbin,fp);
		fwrite(haloprof[i].vsm_this,sizeof(float),nbin,fp);
		fwrite(haloprof[i].vsm_other,sizeof(float),nbin,fp);
		fwrite(haloprof[i].vsm_back,sizeof(float),nbin,fp);
		fwrite(haloprof[i].Um_this,sizeof(float),nbin,fp);
		fwrite(haloprof[i].Um_other,sizeof(float),nbin,fp);
		fwrite(haloprof[i].Um_back,sizeof(float),nbin,fp);
		fwrite(haloprof[i].gvm2_this,sizeof(float),nbin,fp);
		fwrite(haloprof[i].gvm2_other,sizeof(float),nbin,fp);
		fwrite(haloprof[i].gvm2_back,sizeof(float),nbin,fp);
		fwrite(haloprof[i].gvsm_this,sizeof(float),nbin,fp);
		fwrite(haloprof[i].gvsm_other,sizeof(float),nbin,fp);
		fwrite(haloprof[i].gvsm_back,sizeof(float),nbin,fp);
	}
	fclose(fp);
}


/*========to be improved======
 * interpolation or fitting to improve accuracy on radius
 * density calculation of individual particles and then average to improve density accuracy
 * small scale bin to logscale?
 * req to average background inside that radius rather than inside rmax
 * add maximum circular velocity
 * 
 * linkedlist NP_DM and NP_GAS seperately
 * ===*/
