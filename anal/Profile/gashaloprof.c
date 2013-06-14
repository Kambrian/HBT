//use logarithmic radial bin
//the result is no better than linear bin since inner bins have large fluctuation
//perhaps better performance for req?
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

#define NDIV 256
LINKLIST ll;
#ifdef PERIODIC_BDR
#define MAKELL(pos,np,ndiv) make_linklist_box(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_shift_gridid(x,NDIV)
#else
#define MAKELL(pos,np,ndiv) make_linklist(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_round_gridid(x,NDIV)
#endif

GASHALOCAT Cat;
CATALOGUE DMCat;
SUBCATALOGUE SubCat;
typedef struct 
{
	int nbin;
	float rmax;
	float rmin;//min radius of radial log-bins; equals to max(Factor_RMIN*rmax,SofteningHalo); note this is not saved.
	float rvir;//estimated vir for DMhalo
	float req_bk_1;//set to 0 when no bin,to rmax when no back
	float req_all_1;//set to 0 when no bin,to rmax when no back+other
	float req_bk_02;
	float req_all_02;
	float CoM[3];
	float Rvir[3];//DMhalo virial radii,note this is not saved here
	int Mvir_this[3];//halo gas mass inside virial radii
	int Mvir_all[3];//all gas mass inside virial radii
	int mass;
	
	int *n_this;
	int *n_other;
	int *n_back;
	float *r_this;
	float *r_other;
	float *r_back;
} HALOPROFILE;

void halo_radius_eq_global(HALOPROFILE *haloprof);//consider improving accuracy by interpolation or fitting
void halo_bin(int grpid,HALOPROFILE *haloprof);
int fill_haloprof(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
									     int np_this,int np_other,int np_back);
void halo_virial_mass(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
							  int np_this,int np_other,int np_back);										 
int histr_log(float *rdata,int ndata,float rmin,float rmax,int *count,float *rbin,int nbin);
void allocate_haloprof(HALOPROFILE *haloprof,int nbin);
void free_haloprof(HALOPROFILE *haloprof);
void save_gashalo_size(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir);
void save_gashalo_prof(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir);
void load_halo_virial_radii(HALOPROFILE *haloprof,int Ngrps,int Nsnap,char *profdir);

int main(int argc,char **argv)
{
	char outputdir[1024];
	
int Nsnap=99,i;
int grpid,hostsubid,hostnbin,snbin,*sn_host,*n_host,icen;
float dmr,vratio,hostrmax,hostrmin;
double hostspan;
HALOPROFILE *haloprof;

Nsnap=atoi(argv[1]);

logfile=stdout;
sprintf(outputdir,"%s/profile/logbin",SUBCAT_DIR);	

	load_gas_data(Nsnap,SNAPSHOT_DIR);
	load_gashalocat(Nsnap,&Cat,GASCAT_DIR);
	load_group_catalogue(Nsnap,&DMCat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	#ifndef GRPINPUT_INDEX
	fresh_gasID2Index(&Cat,-1); 	//fofcat of JING's data are originally PIndex rather than PID
	#endif
	Cat.ID2Halo=mymalloc(sizeof(int)*NP_GAS);
	prepare_gasind2halo(&Cat);

	MAKELL(Gdat.Pos,NP_GAS,NDIV);
	printf("finished linkedlist\n");

haloprof=mymalloc(sizeof(HALOPROFILE)*Cat.Ngroups);
load_halo_virial_radii(haloprof,Cat.Ngroups,Nsnap,outputdir);
#pragma omp parallel for
for(grpid=0;grpid<Cat.Ngroups;grpid++)
{	
	halo_bin(grpid,haloprof+grpid);//make density profile bins and finding rmax
	halo_radius_eq_global(haloprof+grpid);//calculating req
}
free_linklist(&ll);
save_gashalo_size(haloprof,Cat.Ngroups,Nsnap,outputdir);
save_gashalo_prof(haloprof,Cat.Ngroups,Nsnap,outputdir);
for(grpid=0;grpid<Cat.Ngroups;grpid++)
	free_haloprof(haloprof+grpid);
free(haloprof);
free(Cat.ID2Halo);
free_gashalocat(&Cat);
free_catalogue(&DMCat);
for(i=0;i<SubCat.Nsubs;i++)
	if(SubCat.SubLen[i])
		free(SubCat.PSubArr[i]);
free_sub_catalogue(&SubCat);
return 0;
}
void halo_radius_eq_global(HALOPROFILE *haloprof)//consider improving accuracy by interpolation or fitting
{/*the radius where subdens equals to global  background
* (background averaged inside rmax,not inside rtidal,
* since inside rtidal the background level can be too low)*/
	float *diff_dens_back,*diff_dens_all,av_dens_back,av_dens_all,*vol;
	float rmin,rmax,vratbase;
	double span;
	int i,nbin,sum_back,sum_other,binmax_sub;
	if(nbin=haloprof->nbin)
	{
		sum_back=0;sum_other=0;
		for(i=0;i<nbin;i++)
		{
			sum_back+=haloprof->n_back[i];
			sum_other+=haloprof->n_other[i];
		}
		if(sum_back+sum_other)
		{
			rmin=haloprof->rmin;
			rmax=haloprof->rmax;
			span=log(rmax/rmin);
			vratbase=(exp(3.*span/nbin)-1)/(exp(3.*span)-1);
			vol=mymalloc(sizeof(float)*nbin);
			for(i=0;i<nbin;i++)
				vol[i]=3.*i*span/nbin;
				
			diff_dens_all=mymalloc(sizeof(float)*nbin);	
			av_dens_all=(float)(sum_back+sum_other)*vratbase;
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=logf((float)haloprof->n_this[i]/av_dens_all)-vol[i]; //  logf(/exp(vol[i]))=-vol[i]
			binmax_sub=Fmax_of_vec(diff_dens_all,nbin);
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(diff_dens_all[i]);
			haloprof->req_all_1=haloprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			av_dens_all*=0.02;
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(logf((float)haloprof->n_this[i]/av_dens_all)-vol[i]);
			haloprof->req_all_02=haloprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			free(diff_dens_all);
			
			if(sum_back)
			{
				diff_dens_back=mymalloc(sizeof(float)*nbin);
				av_dens_back=(float)sum_back*vratbase;
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((float)haloprof->n_this[i]/av_dens_back)-vol[i]);
				haloprof->req_bk_1=haloprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
				av_dens_back*=0.02;
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((float)haloprof->n_this[i]/av_dens_back)-vol[i]);
				haloprof->req_bk_02=haloprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
				free(diff_dens_back);
			}
			else
			{
				haloprof->req_bk_1=haloprof->rmax;
				haloprof->req_bk_02=haloprof->rmax;
			}			
			
			free(vol);
		}
		else
		{
			haloprof->req_bk_1=haloprof->rmax;
			haloprof->req_all_1=haloprof->rmax;
			haloprof->req_bk_02=haloprof->rmax;
			haloprof->req_all_02=haloprof->rmax;
		}
	}
	else
	{
		haloprof->req_bk_1=0;
		haloprof->req_all_1=0;
		haloprof->req_bk_02=0;
		haloprof->req_all_02=0;
	}
}

void halo_bin(int grpid,HALOPROFILE *haloprof)
{
	float *cen,rmax,dr;
	float *rp_this,*rp_other,*rp_back;
	int    np_this,np_other,np_back,
		   np_this_max,np_other_max,np_back_max;
	int i,j,k,pid;
	int subbox_grid[3][2];
	
	haloprof->mass=Cat.Len[grpid];
	if(Cat.Len[grpid])
	{
	rp_this=mymalloc(sizeof(float)*Cat.Len[grpid]);
	rp_other=mymalloc(sizeof(float)*Cat.Len[grpid]);
	rp_back=mymalloc(sizeof(float)*Cat.Len[grpid]);
	np_this=np_other=np_back=0;
	np_this_max=np_other_max=np_back_max=Cat.Len[grpid];
	//~ halo_CoM(grpid,haloprof);
	center_of_mass(haloprof->CoM,Cat.PIDorIndex+Cat.Offset[grpid],Cat.Len[grpid],Gdat.Pos);
	if(SubCat.SubLen[SubCat.GrpOffset_Sub[grpid]])
		cen=SubCat.Property[SubCat.GrpOffset_Sub[grpid]].CoM;
	else
		cen=haloprof->CoM;
	haloprof->rvir=comoving_virial_radius(DMCat.Len[grpid]);//use dm mass as halo mass
	rmax=haloprof->rvir*Factor_relax;	
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rmax-ll.range[i][0])/ll.step[i]);
	subbox_grid[i][1]=floor((cen[i]+rmax-ll.range[i][0])/ll.step[i]);
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc(&ll,FIXGRID(i),FIXGRID(j),FIXGRID(k));
				while(pid>=0)
				{
					dr=distance(Gdat.Pos[pid],cen);
					if(dr<rmax)
					{
						if(Cat.ID2Halo[pid]==grpid)
						{
							if(np_this>=np_this_max)
							{
							np_this_max*=2;	
							rp_this=realloc(rp_this,sizeof(float)*np_this_max);
							}
							rp_this[np_this]=dr;
							np_this++;
						}
						else if(Cat.ID2Halo[pid]==-1)
						{
							if(np_back>=np_back_max)
							{
							np_back_max*=2;	
							rp_back=realloc(rp_back,sizeof(float)*np_back_max);
							}
							rp_back[np_back]=dr;
							np_back++;
						}
						else
						{
							if(np_other>=np_other_max)
							{
							np_other_max*=2;	
							rp_other=realloc(rp_other,sizeof(float)*np_other_max);
							}
							rp_other[np_other]=dr;
							np_other++;
						}						
					}
					pid=ll.list[pid];
				}
			}
	fill_haloprof(haloprof,rp_this,rp_other,rp_back,np_this,np_other,np_back);
	halo_virial_mass(haloprof,rp_this,rp_other,rp_back,np_this,np_other,np_back);
	free(rp_this);
	free(rp_other);
	free(rp_back);
	}
	else
	{
	haloprof->CoM[0]=haloprof->CoM[1]=haloprof->CoM[2]=0.;
	haloprof->Mvir_this[0]=haloprof->Mvir_this[1]=haloprof->Mvir_this[2]=0.;
	haloprof->Mvir_all[0]=haloprof->Mvir_all[1]=haloprof->Mvir_all[2]=0.;
	haloprof->nbin=0;
	haloprof->rmax=0;
	haloprof->rmin=0;
	}
}
void halo_virial_mass(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
							  int np_this,int np_other,int np_back)
{
	int i,j,Mvir,Mvir_all;
	float *Rvir;
	Rvir=haloprof->Rvir;

	for(i=0;i<3;i++)
	{
		Mvir=0;
		for(j=0;j<np_this;j++)
			if(rp_this[j]<=Rvir[i])
			Mvir++;
		Mvir_all=Mvir;			
		for(j=0;j<np_other;j++)
			if(rp_other[j]<=Rvir[i])
			Mvir_all++;
		for(j=0;j<np_back;j++)
			if(rp_back[j]<=Rvir[i])
			Mvir_all++;	
	haloprof->Mvir_this[i]=Mvir;
	haloprof->Mvir_all[i]=Mvir_all;
	}
	
}

int fill_haloprof(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
									     int np_this,int np_other,int np_back)
{
	float rmax,rmin;
	int i,nbin,flag_rebin;

	nbin=NBIN_MAX;
	haloprof->n_this=mymalloc(sizeof(int)*nbin);
	haloprof->r_this=mymalloc(sizeof(float)*nbin);
	rmax=rp_this[Fmax_of_vec(rp_this,np_this)];
	haloprof->rmax=rmax;
	rmin=rmax*Factor_RMIN;
	if(rmin<SofteningHalo)
		rmin=SofteningHalo;
	haloprof->rmin=rmin;
	flag_rebin=1;nbin++;
	while(flag_rebin&&nbin>NBIN_MIN)
	{
	nbin--;
	flag_rebin=histr_log(rp_this,np_this,rmin,rmax,haloprof->n_this,haloprof->r_this,nbin);
	}
	if(NBIN_MIN==nbin)
	{
	free(haloprof->n_this);haloprof->n_this=NULL;
	free(haloprof->r_this);haloprof->r_this=NULL;
	haloprof->nbin=0;
	return 0;
	}
	haloprof->nbin=nbin;
	haloprof->n_this=realloc(haloprof->n_this,sizeof(int)*nbin);
	haloprof->r_this=realloc(haloprof->r_this,sizeof(float)*nbin);
	haloprof->n_other=mymalloc(sizeof(int)*nbin);
	haloprof->r_other=mymalloc(sizeof(float)*nbin);
	haloprof->n_back=mymalloc(sizeof(int)*nbin);
	haloprof->r_back=mymalloc(sizeof(float)*nbin);
	histr_log(rp_other,np_other,rmin,rmax,haloprof->n_other,haloprof->r_other,nbin);
	histr_log(rp_back,np_back,rmin,rmax,haloprof->n_back,haloprof->r_back,nbin);
	return nbin;
}

int histr_log(float *rdata,int ndata,float rmin,float rmax,int *count,float *rbin,int nbin)
{
	//return 1 if have empty bins,0 otherwise
	int i,dr_bin,flag_emptybin;
	double span_order;
	for(i=0;i<nbin;i++)
	{
		count[i]=0;
		rbin[i]=0;
	}
	span_order=log(rmax/rmin);
	for(i=0;i<ndata;i++)
	{
		if(rdata[i]>=rmin&&rdata[i]<rmax)
		{
		dr_bin=floor(log(rdata[i]/rmin)/span_order*nbin);
		if(dr_bin<0) dr_bin=0;
		if(dr_bin>=nbin) dr_bin=nbin-1;
		count[dr_bin]++;
		rbin[dr_bin]+=rdata[i];
		}
	}
	flag_emptybin=0;
	for(i=0;i<nbin;i++)
	{
		if(count[i])
		rbin[i]/=count[i];
		else
		{
		flag_emptybin=1;
		rbin[i]=(exp((i+1)*span_order/nbin)+exp(i*span_order/nbin))*0.5*rmin;
		}
	}
	return flag_emptybin;	
}
void allocate_haloprof(HALOPROFILE *haloprof,int nbin)
{
	haloprof->nbin=nbin;
	if(nbin)
	{
		haloprof->n_this=mymalloc(sizeof(int)*nbin);
		haloprof->n_other=mymalloc(sizeof(int)*nbin);
		haloprof->n_back=mymalloc(sizeof(int)*nbin);
		haloprof->r_this=mymalloc(sizeof(float)*nbin);
		haloprof->r_other=mymalloc(sizeof(float)*nbin);
		haloprof->r_back=mymalloc(sizeof(float)*nbin);
	}
}

void free_haloprof(HALOPROFILE *haloprof)
{
	if(haloprof->nbin)
	{
		free(haloprof->n_this);
		free(haloprof->n_other);
		free(haloprof->n_back);
		free(haloprof->r_this);
		free(haloprof->r_other);
		free(haloprof->r_back);
	}
}

void save_gashalo_size(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/ghalo_size_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Nsubs;i++)
	{
		fwrite(&haloprof[i].nbin,sizeof(int),1,fp);
		fwrite(&haloprof[i].rmax,sizeof(float),1,fp);
		fwrite(&haloprof[i].rvir,sizeof(float),1,fp);
		fwrite(&haloprof[i].req_bk_1,sizeof(float),1,fp);
		fwrite(&haloprof[i].req_all_1,sizeof(float),1,fp);
		fwrite(&haloprof[i].req_bk_02,sizeof(float),1,fp);
		fwrite(&haloprof[i].req_all_02,sizeof(float),1,fp);
		fwrite(haloprof[i].CoM,sizeof(float),3,fp);
		fwrite(haloprof[i].Mvir_this,sizeof(int),3,fp);
		fwrite(haloprof[i].Mvir_all,sizeof(int),3,fp);
		fwrite(&haloprof[i].mass,sizeof(int),1,fp);
	}
	fclose(fp);
}
void save_gashalo_prof(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i,nbin;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/ghalo_prof_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Nsubs;i++)
	{
		nbin=haloprof[i].nbin;
		fwrite(haloprof[i].n_this,sizeof(int),nbin,fp);
		fwrite(haloprof[i].n_other,sizeof(int),nbin,fp);
		fwrite(haloprof[i].n_back,sizeof(int),nbin,fp);
		fwrite(haloprof[i].r_this,sizeof(float),nbin,fp);
		fwrite(haloprof[i].r_other,sizeof(float),nbin,fp);
		fwrite(haloprof[i].r_back,sizeof(float),nbin,fp);
	}
	fclose(fp);
}

void load_halo_virial_radii(HALOPROFILE *haloprof,int Ngrps,int Nsnap,char *profdir)
{
	int i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/halo_size_%03d",profdir,Nsnap);
	myfopen(fp,buf,"r");
	fseek(fp,17*4L,SEEK_SET);
	fread(haloprof[0].Rvir,sizeof(float),3,fp);
	for(i=1;i<Ngrps;i++)
	{
		fseek(fp,21*4L,SEEK_CUR);
		fread(haloprof[i].Rvir,sizeof(float),3,fp);
	}
	fclose(fp);
}
/*========to be improved======
 * interpolation or fitting to improve accuracy on radius
 * density calculation of individual particles and then average to improve density accuracy
 * small scale bin to logscale?
 * req to average background inside that radius rather than inside rmax
 * add maximum circular velocity
 * ===*/
