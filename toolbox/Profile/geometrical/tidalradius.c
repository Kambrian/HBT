#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define NBIN_MAX 200
#define NBIN_MIN 5
#define Factor_relax 2

#define NDIV 200
LINKLIST ll;
#ifdef PERIODIC_BDR
#define MAKELL(pos,np,ndiv) make_linklist_box(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_shift_gridid(x,NDIV)
#else
#define MAKELL(pos,np,ndiv) make_linklist(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_round_gridid(x,NDIV)
#endif

//to use most bound center
#define SUBCOM(i) Pdat.Pos[SubCat.PSubArr[i][0]]  
//to use core com
//~ #define SUBCOM(i) SubCat.Property[i].CoM

//~ #define SUBFIND_DIR "/home/kambrain/data/6702DM/subcatS"   //define this properly to apply on subfind, 
															//undef for HBT

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

int *PartHost;
SUBCATALOGUE SubCat;

typedef struct 
{
	int nbin;
	int flag_badbins;//1 means badbin,-1 means bad host,0 normal,2 means rcen>rmax
	float rmax;
	float rvir;
	float req_bk_1;//set to 0 when no bin,to rmax when no back
	float req_all_1;//set to 0 when no bin,to rmax when no back+other
	float req_bk_02;
	float req_all_02;
	float rtidal;//set to rvir when cannot determine (e.g.,mainsubs/quasi-subs and subhalos inside a bad-host)
	float rcen;
	
	int *n_this;
	int *n_other;
	int *n_back;
	float *r_this;
	float *r_other;
	float *r_back;
} SUBPROFILE;

void sub_radius_eq_local(SUBPROFILE *subprof);//consider improving accuracy by interpolation or fitting
void sub_radius_eq_global(SUBPROFILE *subprof);//consider improving accuracy by interpolation or fitting
void sub_bin(int subid,SUBPROFILE *subprof);
int fill_subprof(SUBPROFILE *subprof,float *rp_this,float *rp_other,float *rp_back,
									     int np_this,int np_other,int np_back);
int histr(float *rdata,int ndata,float rmax,int *count,float *rbin,int nbin);
void allocate_subprof(SUBPROFILE *subprof,int nbin);
void free_subprof(SUBPROFILE *subprof);
void save_subs_size(SUBPROFILE *subprof,int Nsubs,int Nsnap,char *outputdir);
void save_subs_prof(SUBPROFILE *subprof,int Nsubs,int Nsnap,char *outputdir);

int main()
{
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	
int Nsnap=99,i;
int subid,grpid,hostsubid,hostnbin,snbin,*sn_host,*n_host,icen;
float dmr,vratio,hostrmax;
SUBPROFILE *subprof;

logfile=stdout;
sprintf(outputdir,"%s/profile",SUBCAT_DIR);	

load_sub_catalogue(Nsnap,&SubCat,inputdir);
load_particle_data(Nsnap,snapdir);
fresh_ID2Index(&SubCat,-2);

PartHost=mymalloc(sizeof(int)*NP_DM);
for(i=0;i<NP_DM;i++)
	PartHost[i]=-1;
for(subid=0;subid<SubCat.Nsubs;subid++)
{
	for(i=0;i<SubCat.SubLen[subid];i++)
		PartHost[SubCat.PSubArr[subid][i]]=subid;
}

	MAKELL(Pdat.Pos,NP_DM,NDIV);
	printf("finished linkedlist\n");

subprof=mymalloc(sizeof(SUBPROFILE)*SubCat.Nsubs);
for(subid=0;subid<SubCat.Nsubs;subid++)
{	
	sub_bin(subid,subprof+subid);//make density profile bins and finding rmax
	hostsubid=(SubCat.HaloChains[subid].HostID<0)?(subid):(SubCat.GrpOffset_Sub[SubCat.HaloChains[subid].HostID]);
	subprof[subid].rcen=(SubCat.SubLen[subid]>0)?distance(SUBCOM(subid),SUBCOM(hostsubid)):0;//finding rcen
	sub_radius_eq_global(subprof+subid);//calculating req
	subprof[subid].flag_badbins=0;
}
	free_linklist(&ll);
//calculating tidal radius
for(grpid=0;grpid<SubCat.Ngroups;grpid++)
{
	hostsubid=SubCat.GrpOffset_Sub[grpid];
	subprof[hostsubid].rtidal=subprof[hostsubid].rvir;
	if(SubCat.GrpLen_Sub[grpid]>1)
	{
		if(hostnbin=subprof[hostsubid].nbin)
		{
			hostrmax=subprof[hostsubid].rmax;
			//~ n_host=subprof[hostsubid].n_this;
			n_host=mymalloc(sizeof(int)*hostnbin);//all as tidal
			sn_host=mymalloc(sizeof(int)*hostnbin);
			snbin=0;
			for(i=0;i<hostnbin;i++)
			{
				n_host[i]=subprof[hostsubid].n_this[i]+subprof[hostsubid].n_other[i]+subprof[hostsubid].n_back[i];//all as tidal
				snbin+=n_host[i];
				sn_host[i]=snbin;
			}
			for(subid=hostsubid+1;subid<hostsubid+SubCat.GrpLen_Sub[grpid];subid++)
			{
				icen=floor(subprof[subid].rcen/hostrmax*hostnbin);
				if(icen>=hostnbin)//a far-away sub
				{	
				subprof[subid].flag_badbins=2;//more than outer	
				icen=hostnbin-1;
				dmr=0.;
				}
				else
				{
				vratio=(float)icen/(icen+1);
				vratio=vratio*vratio*vratio;
				dmr=3.*n_host[icen]/(1-vratio)/sn_host[icen];
				}
				if(dmr>=3)//probably a badly sampled or a quite loose host (badbin)
				{
				subprof[subid].flag_badbins=1;
				subprof[subid].rtidal=subprof[subid].rvir;
				}
				else
				subprof[subid].rtidal=subprof[subid].rcen*pow(SubCat.SubLen[subid]/(3.-dmr)/sn_host[icen],1.0/3);
			}
			free(n_host);//all as tidal
			free(sn_host);
		}
		else//bad_host
		{
			subprof[subid].flag_badbins=-1;
			subprof[subid].rtidal=subprof[subid].rvir;
		}
	}
}
for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++)//quasi
	subprof[subid].rtidal=subprof[subid].rvir;
save_subs_size(subprof,SubCat.Nsubs,Nsnap,outputdir);
save_subs_prof(subprof,SubCat.Nsubs,Nsnap,outputdir);
for(subid=0;subid<SubCat.Nsubs;subid++)
	free_subprof(subprof+subid);
free(subprof);
free(PartHost);	
free_sub_catalogue(&SubCat);

return 0;
}
void sub_radius_eq_global(SUBPROFILE *subprof)//consider improving accuracy by interpolation or fitting
{/*the radius where subdens equals to global  background
* (background averaged inside rmax,not inside rtidal,
* since inside rtidal the background level can be too low)*/
	float *diff_dens_back,*diff_dens_all,av_dens_back,av_dens_all;
	int i,nbin,*vol,sum_back,sum_other,binmax_sub;
	if(nbin=subprof->nbin)
	{
		sum_back=0;sum_other=0;
		for(i=0;i<nbin;i++)
		{
			sum_back+=subprof->n_back[i];
			sum_other+=subprof->n_other[i];
		}
		if(sum_back+sum_other)
		{
			vol=mymalloc(sizeof(int)*nbin);
			for(i=1;i<=nbin;i++)
				vol[i-1]=i*i*i;
			for(i=1;i<nbin;i++)
				vol[i]=vol[i]-vol[i-1];
				
			diff_dens_all=mymalloc(sizeof(float)*nbin);	
			av_dens_all=(float)(sum_back+sum_other)/(nbin*nbin*nbin);
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=logf((float)subprof->n_this[i]/vol[i]/av_dens_all);
			binmax_sub=Fmax_of_vec(diff_dens_all,nbin);
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(diff_dens_all[i]);
			subprof->req_all_1=subprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			av_dens_all*=0.02;
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(logf((float)subprof->n_this[i]/vol[i]/av_dens_all));
			subprof->req_all_02=subprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			free(diff_dens_all);
			
			if(sum_back)
			{
				diff_dens_back=mymalloc(sizeof(float)*nbin);
				av_dens_back=(float)sum_back/(nbin*nbin*nbin);
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((float)subprof->n_this[i]/vol[i]/av_dens_back));
				subprof->req_bk_1=subprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
				av_dens_back*=0.02;
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((float)subprof->n_this[i]/vol[i]/av_dens_back));
				subprof->req_bk_02=subprof->r_this[binmax_sub+Fmin_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
				free(diff_dens_back);
			}
			else
			{
				subprof->req_bk_1=subprof->rmax;
				subprof->req_bk_02=subprof->rmax;
			}			
			
			free(vol);
		}
		else
		{
			subprof->req_bk_1=subprof->rmax;
			subprof->req_all_1=subprof->rmax;
			subprof->req_bk_02=subprof->rmax;
			subprof->req_all_02=subprof->rmax;
		}
	}
	else
	{
		subprof->req_bk_1=0;
		subprof->req_all_1=0;
		subprof->req_bk_02=0;
		subprof->req_all_02=0;
	}
}
void sub_radius_eq_local(SUBPROFILE *subprof)//the radius where subdens equal to local background
{
	float *diff_dens_back,*diff_dens_all;
	int i;
	if(subprof->nbin)
	{
	diff_dens_back=mymalloc(sizeof(float)*subprof->nbin);
	diff_dens_all=mymalloc(sizeof(float)*subprof->nbin);

	for(i=0;i<subprof->nbin;i++)
	{
		diff_dens_back[i]=fabs(subprof->n_this[i]-subprof->n_back[i])/(float)subprof->n_this[i];
		diff_dens_all[i]=fabs(subprof->n_this[i]-subprof->n_back[i]-subprof->n_other[i])/(float)subprof->n_this[i];
	}
	subprof->req_bk_1=subprof->r_this[Fmin_of_vec(diff_dens_back,subprof->nbin)];
	subprof->req_all_1=subprof->r_this[Fmin_of_vec(diff_dens_all,subprof->nbin)];
	
	for(i=0;i<subprof->nbin;i++)
	{
		diff_dens_back[i]=fabs(subprof->n_this[i]-0.02*subprof->n_back[i])/(float)subprof->n_this[i];
		diff_dens_all[i]=fabs(subprof->n_this[i]-0.02*(subprof->n_back[i]+subprof->n_other[i]))/(float)subprof->n_this[i];
	}
	subprof->req_bk_02=subprof->r_this[Fmin_of_vec(diff_dens_back,subprof->nbin)];
	subprof->req_all_02=subprof->r_this[Fmin_of_vec(diff_dens_all,subprof->nbin)];
	free(diff_dens_back);
	free(diff_dens_all);
	}
	else
	{
	subprof->req_bk_1=0;
	subprof->req_all_1=0;
	subprof->req_bk_02=0;
	subprof->req_all_02=0;
	}
}

void sub_bin(int subid,SUBPROFILE *subprof)
{
	float *cen,rmax,dr;
	float *rp_this,*rp_other,*rp_back;
	int    np_this,np_other,np_back,
		   np_this_max,np_other_max,np_back_max;
	int i,j,k,pid;
	int subbox_grid[3][2];
	
	if(SubCat.SubLen[subid])
	{
	rp_this=mymalloc(sizeof(float)*SubCat.SubLen[subid]);
	rp_other=mymalloc(sizeof(float)*SubCat.SubLen[subid]);
	rp_back=mymalloc(sizeof(float)*SubCat.SubLen[subid]);
	np_this=np_other=np_back=0;
	np_this_max=np_other_max=np_back_max=SubCat.SubLen[subid];
	//cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
	cen=SUBCOM(subid);
	subprof->rvir=comoving_virial_radius(SubCat.SubLen[subid]);
	rmax=subprof->rvir*Factor_relax;	
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
						if(PartHost[pid]==subid)
						{
							if(np_this>=np_this_max)
							{
							np_this_max*=2;	
							rp_this=realloc(rp_this,sizeof(float)*np_this_max);
							}
							rp_this[np_this]=dr;
							np_this++;
						}
						else if(PartHost[pid]==-1)
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
	fill_subprof(subprof,rp_this,rp_other,rp_back,np_this,np_other,np_back);
	free(rp_this);
	free(rp_other);
	free(rp_back);
	}
	else
	{
	subprof->nbin=0;
	subprof->rmax=0;
	}
}
int fill_subprof(SUBPROFILE *subprof,float *rp_this,float *rp_other,float *rp_back,
									     int np_this,int np_other,int np_back)
{
	float rmax;
	int i,nbin,flag_rebin;

	nbin=NBIN_MAX;
	subprof->n_this=mymalloc(sizeof(int)*nbin);
	subprof->r_this=mymalloc(sizeof(float)*nbin);
	rmax=rp_this[Fmax_of_vec(rp_this,np_this)];
	subprof->rmax=rmax;
	flag_rebin=1;nbin++;
	while(flag_rebin&&nbin>NBIN_MIN)
	{
	nbin--;
	flag_rebin=histr(rp_this,np_this,rmax,subprof->n_this,subprof->r_this,nbin);
	}
	if(NBIN_MIN==nbin)
	{
	free(subprof->n_this);subprof->n_this=NULL;
	free(subprof->r_this);subprof->r_this=NULL;
	subprof->nbin=0;
	return 0;
	}
	subprof->nbin=nbin;
	subprof->n_this=realloc(subprof->n_this,sizeof(int)*nbin);
	subprof->r_this=realloc(subprof->r_this,sizeof(float)*nbin);
	subprof->n_other=mymalloc(sizeof(int)*nbin);
	subprof->r_other=mymalloc(sizeof(float)*nbin);
	subprof->n_back=mymalloc(sizeof(int)*nbin);
	subprof->r_back=mymalloc(sizeof(float)*nbin);
	histr(rp_other,np_other,rmax,subprof->n_other,subprof->r_other,nbin);
	histr(rp_back,np_back,rmax,subprof->n_back,subprof->r_back,nbin);
	return nbin;
}

int histr(float *rdata,int ndata,float rmax,int *count,float *rbin,int nbin)
{
	//return 1 if have empty bins,0 otherwise
	int i,dr_bin,flag_emptybin;
	for(i=0;i<nbin;i++)
	{
		count[i]=0;
		rbin[i]=0;
	}
	for(i=0;i<ndata;i++)
	{
		dr_bin=floor(rdata[i]/rmax*nbin);
		if(dr_bin<nbin)
		{
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
		rbin[i]=(i+0.5)/nbin*rmax;
		}
	}
	return flag_emptybin;	
}
void allocate_subprof(SUBPROFILE *subprof,int nbin)
{
	subprof->nbin=nbin;
	if(nbin)
	{
		subprof->n_this=mymalloc(sizeof(int)*nbin);
		subprof->n_other=mymalloc(sizeof(int)*nbin);
		subprof->n_back=mymalloc(sizeof(int)*nbin);
		subprof->r_this=mymalloc(sizeof(float)*nbin);
		subprof->r_other=mymalloc(sizeof(float)*nbin);
		subprof->r_back=mymalloc(sizeof(float)*nbin);
	}
}

void free_subprof(SUBPROFILE *subprof)
{
	if(subprof->nbin)
	{
		free(subprof->n_this);
		free(subprof->n_other);
		free(subprof->n_back);
		free(subprof->r_this);
		free(subprof->r_other);
		free(subprof->r_back);
	}
}

void save_subs_size(SUBPROFILE *subprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/sub_size_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Nsubs;i++)
	{
		fwrite(&subprof[i].nbin,sizeof(int),1,fp);
		fwrite(&subprof[i].flag_badbins,sizeof(int),1,fp);
		fwrite(&subprof[i].rmax,sizeof(float),1,fp);
		fwrite(&subprof[i].rvir,sizeof(float),1,fp);
		fwrite(&subprof[i].req_bk_1,sizeof(float),1,fp);
		fwrite(&subprof[i].req_all_1,sizeof(float),1,fp);
		fwrite(&subprof[i].req_bk_02,sizeof(float),1,fp);
		fwrite(&subprof[i].req_all_02,sizeof(float),1,fp);
		fwrite(&subprof[i].rtidal,sizeof(float),1,fp);
		fwrite(&subprof[i].rcen,sizeof(float),1,fp);
	}
	fclose(fp);
}
void save_subs_prof(SUBPROFILE *subprof,int Nsubs,int Nsnap,char *outputdir)
{
	int i,nbin;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/sub_prof_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Nsubs;i++)
	{
		nbin=subprof[i].nbin;
		fwrite(subprof[i].n_this,sizeof(int),nbin,fp);
		fwrite(subprof[i].n_other,sizeof(int),nbin,fp);
		fwrite(subprof[i].n_back,sizeof(int),nbin,fp);
		fwrite(subprof[i].r_this,sizeof(float),nbin,fp);
		fwrite(subprof[i].r_other,sizeof(float),nbin,fp);
		fwrite(subprof[i].r_back,sizeof(float),nbin,fp);
	}
	fclose(fp);
}


/*========to be improved======
 * interpolation or fitting to improve accuracy on radius
 * density calculation of individual particles and then average to improve density accuracy
 * small scale bin to logscale?
 * req to average background inside that radius rather than inside rmax
 * ===*/
