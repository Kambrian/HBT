//use logarithmic radial bin
//the result is no better than linear bin since inner bins have large fluctuation
//perhaps better performance for req?
//28-11-2013: extended to run without subcat.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define NBIN_FIX 20
#define VIR_TH 0
#define VIR_C200 1
#define VIR_B200 2
#define RBIN_TYPE VIR_C200

#define NBIN_MAX 60
#define NBIN_MIN 5
#define Factor_RMIN 1e-2
#define Factor_relax 3

#define NDIV 256
LINKLIST ll;
#ifdef PERIODIC_BDR
#define MAKELL() make_linklist(&ll,NP_DM,NDIV,Pdat.Pos,GetArrPos,1)
#define FIXGRID(x) linklist_shift_gridid(x,NDIV)
#else
#define MAKELL() make_linklist(&ll,NP_DM,NDIV,Pdat.Pos,GetArrPos,0)
#define FIXGRID(x) linklist_round_gridid(x,NDIV)
#endif
CATALOGUE Cat;
#ifndef WITHOUT_SUBCAT //to run without subcat
SUBCATALOGUE SubCat;
//~ #define SUBFIND_DIR "/home/kambrain/data/6402U100/subcatS"
#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif
#endif

typedef struct 
{
	int nbin;
	float rmax;
	float rmin;//min radius of radial log-bins; equals to max(Factor_RMIN*rmax,SofteningHalo); note this is not saved.
	float rvir;
	float req_bk_1;//set to 0 when no bin,to rmax when no back
	float req_all_1;//set to 0 when no bin,to rmax when no back+other
	float req_bk_02;
	float req_all_02;
	float CoM[3];//halo CoM,for comparison
	float Cen[3];//Center of mainsub,the profiles are w.r.t this center
	int mass;//fof mass
	int Mvir[3];
	float Rvir[3];//[tophat,c200,b200],comoving
	int flag_badvir[3];
	int flag_fakehalo;//set to 1 when halo is not self-bound,to 0 otherwise.
	
	int *n_this;
	int *n_other;
	int *n_back;
	float *r_this;
	float *r_other;
	float *r_back;
} HALOPROFILE;

void halo_radius_eq_global(HALOPROFILE *haloprof);//consider improving accuracy by interpolation or fitting
void halo_bin(int grpid,HALOPROFILE *haloprof,float virialF[3]);
void halo_virial_factor(float virialF[3]);
void halo_virial_radius(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
							  int np_this,int np_other,int np_back,float virialF[3]);
int fill_haloprof(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
									     int np_this,int np_other,int np_back);
int histr_log(float *rdata,int ndata,float rmin,float rmax,int *count,float *rbin,int nbin);
void allocate_haloprof(HALOPROFILE *haloprof,int nbin);
void free_haloprof(HALOPROFILE *haloprof);
void save_halos_size(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir);
void save_halos_prof(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir);

int main(int argc,char **argv)
{
char outputdir[1024];
	
int Nsnap=99,i;
int grpid;
float virialF[3];
HALOPROFILE *haloprof;

	if(argc!=2)
	{
		printf("usage: %s [Nsnap]\n",argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);

	logfile=stdout;
	sprintf(outputdir,"%s/profile",SUBCAT_DIR);	
	mkdir(outputdir,0755);
	sprintf(outputdir,"%s/profile/logbin",SUBCAT_DIR);	
	mkdir(outputdir,0755);

	load_particle_data(Nsnap,SNAPSHOT_DIR);	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
#ifndef WITHOUT_SUBCAT
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
#endif
	fill_PIDHash();
	fresh_ID2Index(&Cat,-1); 	//fofcat of JING's data are originally PIndex rather than PID
	free_PIDHash();
	prepare_ind2halo(&Cat);

	MAKELL();
	printf("finished linkedlist\n");

	haloprof=mymalloc(sizeof(HALOPROFILE)*Cat.Ngroups);
	halo_virial_factor(virialF);
	#pragma omp parallel for schedule(dynamic,1)
	for(grpid=0;grpid<Cat.Ngroups;grpid++)
	{	
		halo_bin(grpid,haloprof+grpid,virialF);//make density profile bins and finding rmax
		halo_radius_eq_global(haloprof+grpid);//calculating req
	}
	
	free_linklist(&ll);
	save_halos_size(haloprof,Cat.Ngroups,Nsnap,outputdir);
	save_halos_prof(haloprof,Cat.Ngroups,Nsnap,outputdir);
	for(grpid=0;grpid<Cat.Ngroups;grpid++)
		free_haloprof(haloprof+grpid);
	free(haloprof);

	free_catalogue(&Cat);
#ifndef WITHOUT_SUBCAT
	#ifndef SUBFIND_DIR
	erase_sub_catalogue(&SubCat);
	#endif
#endif
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

void halo_bin(int grpid,HALOPROFILE *haloprof,float virialF[3])
{
	HBTReal cen[3],rmax,dr;
	float *rp_this,*rp_other,*rp_back;
	int    np_this,np_other,np_back,
		   np_this_max,np_other_max,np_back_max;
	HBTInt i,j,k,pid;
	int subbox_grid[3][2];
	
	haloprof->mass=Cat.Len[grpid];
	//~ if(Cat.Len[grpid])
	//~ {
	rp_this=mymalloc(sizeof(float)*Cat.Len[grpid]);
	rp_other=mymalloc(sizeof(float)*Cat.Len[grpid]);
	rp_back=mymalloc(sizeof(float)*Cat.Len[grpid]);
	np_this=np_other=np_back=0;
	np_this_max=np_other_max=np_back_max=Cat.Len[grpid];
	//cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
	center_of_mass(cen,Cat.PIDorIndex+Cat.Offset[grpid],Cat.Len[grpid],Pdat.Pos);
	for(i=0;i<3;i++) haloprof->CoM[i]=cen[i];
	//~ cen=haloprof->CoM;
#ifndef WITHOUT_SUBCAT
	for(i=0;i<3;i++)
	haloprof->Cen[i]=SubCat.Property[SubCat.GrpOffset_Sub[grpid]].CoM[i];
	if(SubCat.SubLen[SubCat.GrpOffset_Sub[grpid]])
	{
		for(i=0;i<3;i++)
		cen[i]=haloprof->Cen[i];
		haloprof->flag_fakehalo=0;
	}
	else
	{
		for(i=0;i<3;i++)
		cen[i]=haloprof->CoM[i];
		haloprof->flag_fakehalo=1;
	}
#else
	moving_center(Cat.PIDorIndex+Cat.Offset[grpid],Cat.Len[grpid],cen);
	for(i=0;i<3;i++)
	  haloprof->Cen[i]=cen[i];
	haloprof->flag_fakehalo=1;//not test, always set to 1
#endif
	haloprof->rvir=comoving_virial_radius(Cat.Len[grpid]);
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
					dr=distance(Pdat.Pos[pid],cen);
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
	halo_virial_radius(haloprof,rp_this,rp_other,rp_back,np_this,np_other,np_back,virialF);		
	fill_haloprof(haloprof,rp_this,rp_other,rp_back,np_this,np_other,np_back);
	free(rp_this);
	free(rp_other);
	free(rp_back);
	//~ }
	//~ else
	//~ {
	//~ haloprof->CoM[0]=haloprof->CoM[1]=haloprof->CoM[2]=0.;
	//~ haloprof->Cen[0]=haloprof->Cen[1]=haloprof->Cen[2]=0.;
	//~ haloprof->Rvir[0]=haloprof->Rvir[1]=haloprof->Rvir[2]=0.;
	//~ haloprof->Mvir[0]=haloprof->Mvir[1]=haloprof->Mvir[2]=0.;
	//~ haloprof->flag_badvir[0]=haloprof->flag_badvir[1]=haloprof->flag_badvir[2]=0.;
	//~ haloprof->nbin=0;
	//~ haloprof->rmax=0;
	//~ haloprof->rmin=0;
	//~ }
}
void halo_virial_factor(float virialF[3])
{
	float Hratio,scaleF,virialF_tophat,virialF_c200,virialF_b200,x,OmegaZ;
	scaleF=header.time;
	Hratio=header.Hz/HUBBLE0;
	#ifdef OMEGA0
	OmegaZ=OMEGA0/(scaleF*scaleF*scaleF)/Hratio/Hratio;
	#else
	OmegaZ=header.Omega0/(scaleF*scaleF*scaleF)/Hratio/Hratio;
	#endif
	x=OmegaZ-1;
	virialF_tophat=18.0*3.1416*3.1416+82.0*x-39.0*x*x;//<Rho_vir>/Rho_cri
	virialF_c200=200.;
	virialF_b200=200.*OmegaZ;//virialF w.r.t contemporary critical density 
	virialF[0]=virialF_tophat;
	virialF[1]=virialF_c200;
	virialF[2]=virialF_b200;
}
void halo_virial_radius(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
							  int np_this,int np_other,int np_back,float virialF[3])
{
	float tol=1e-5;
	int i,np_new,np_in,np_out;
	float *rp_new,*rp_in,*rp_out,*rp_tmp,rvir,rdiv;
	int virtype;
	
	for(virtype=0;virtype<3;virtype++)
	{
	np_in=np_this+np_other+np_back;
	np_out=0;
	rp_in=mymalloc(sizeof(float)*np_in);
	rp_out=mymalloc(sizeof(float)*np_in);
	rp_new=mymalloc(sizeof(float)*np_in);
	memcpy(rp_in,rp_this,sizeof(float)*np_this);
	memcpy(rp_in+np_this,rp_other,sizeof(float)*np_other);
	memcpy(rp_in+np_this+np_other,rp_back,sizeof(float)*np_back);	
	rvir=pow(2.0*G*np_in*header.mass[1]/virialF[virtype]/header.Hz/header.Hz,1.0/3)/header.time;
	rdiv=rvir*10;//just to make rvir<rdiv to start the loop
	while(1)
	{
		if(rvir<rdiv)
		{
		np_new=0;
		for(i=0;i<np_in;i++)
		{
			if(rp_in[i]>rvir)
			{
			rp_out[np_out]=rp_in[i];
			np_out++;
			}
			else
			{
			rp_new[np_new]=rp_in[i];
			np_new++;
			}
		}
		rp_tmp=rp_new;rp_new=rp_in;rp_in=rp_tmp;
		np_in=np_new;
		}
		else if(rvir>rdiv)
		{
		np_new=0;
		for(i=0;i<np_out;i++)
		{
			if(rp_out[i]<=rvir)
			{
			rp_in[np_in]=rp_out[i];
			np_in++;
			}
			else
			{
			rp_new[np_new]=rp_out[i];
			np_new++;
			}
		}
		rp_tmp=rp_new;rp_new=rp_out;rp_out=rp_tmp;
		np_out=np_new;	
		}
		
		rdiv=rvir;
		rvir=pow(2.0*G*np_in*header.mass[1]/virialF[virtype]/header.Hz/header.Hz,1.0/3)/header.time;
		if(0==np_in)
		{
			haloprof->flag_badvir[virtype]=1;
			break;
		}
		if(0==np_out)
		{
			haloprof->flag_badvir[virtype]=-1;
			break;
		}
		if(fabs((rvir-rdiv)/rvir)<tol)
		{
			haloprof->flag_badvir[virtype]=0;
			break;
		}
	}
	free(rp_in);
	free(rp_out);
	free(rp_new);
	if(rvir>haloprof->rvir*Factor_relax)haloprof->flag_badvir[virtype]=-1;
	haloprof->Rvir[virtype]=rvir;
	haloprof->Mvir[virtype]=np_in;
	}
}


int fill_haloprof(HALOPROFILE *haloprof,float *rp_this,float *rp_other,float *rp_back,
									     int np_this,int np_other,int np_back)
{
	float rmax,rmin;
	int nbin,flag_rebin;

	nbin=NBIN_FIX;
// 	nbin=NBIN_MAX;
	haloprof->n_this=mymalloc(sizeof(int)*nbin);
	haloprof->r_this=mymalloc(sizeof(float)*nbin);
	rmax=haloprof->Rvir[RBIN_TYPE];
// 	rmax=rp_this[Fmax_of_vec(rp_this,np_this)];
	haloprof->rmax=rmax;
// 	rmin=rmax*Factor_RMIN;
// 	if(rmin<SofteningHalo)
		rmin=SofteningHalo;
	haloprof->rmin=rmin;
	histr_log(rp_this,np_this,rmin,rmax,haloprof->n_this,haloprof->r_this,nbin);
// 	flag_rebin=1;nbin++;
// 	while(flag_rebin&&nbin>NBIN_MIN)
// 	{
// 	nbin--;
// 	flag_rebin=histr_log(rp_this,np_this,rmin,rmax,haloprof->n_this,haloprof->r_this,nbin);
// 	}
// 	if(NBIN_MIN==nbin)
// 	{
// 	free(haloprof->n_this);haloprof->n_this=NULL;
// 	free(haloprof->r_this);haloprof->r_this=NULL;
// 	haloprof->nbin=0;
// 	return 0;
// 	}
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

void save_halos_size(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/halo_size_%03d",outputdir,Nsnap);
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
		fwrite(haloprof[i].Cen,sizeof(float),3,fp);
		fwrite(&haloprof[i].mass,sizeof(int),1,fp);
		fwrite(haloprof[i].Mvir,sizeof(int),3,fp);
		fwrite(haloprof[i].Rvir,sizeof(float),3,fp);
		fwrite(haloprof[i].flag_badvir,sizeof(int),3,fp);
		fwrite(&haloprof[i].flag_fakehalo,sizeof(int),1,fp);
	}
	fclose(fp);
}
void save_halos_prof(HALOPROFILE *haloprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i,nbin;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/halo_prof_%03d",outputdir,Nsnap);
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

/*========to be improved======
 * interpolation or fitting to improve accuracy on radius
 * density calculation of individual particles and then average to improve density accuracy
 * small scale bin to logscale?
 * req to average background inside that radius rather than inside rmax
 * add maximum circular velocity
 * ===*/
