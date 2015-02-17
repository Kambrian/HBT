//use logarithmic radial bin
//the result is no better than linear bin since inner bins have large fluctuation
//perhaps better performance for req?
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

//#define AQUA_SIZE

#define TIDAL_MAX 3   //2 for no centrifugal force, 3 for circular orbit

#define GRPLEN_MIN 0
#define SUBLEN_MIN 0
#define GRPID_MAX SubCat.Ngroups
//~ #define GRPID_MAX 100

#define NBIN_MAX 80
#define NBIN_MIN 5
#define Factor_RMIN 1e-2
#define Factor_relax 2

#define NDIV 200
LINKLIST ll;
#ifdef PERIODIC_BDR
#define MAKELL(pos,np,ndiv) make_linklist(&ll,np,ndiv,pos, GetArrPos,1)
#define FIXGRID(x) linklist_shift_gridid(x,NDIV)
#else
#define MAKELL(pos,np,ndiv) make_linklist(&ll,np,ndiv,pos, GetArrPos,0)
#define FIXGRID(x) linklist_round_gridid(x,NDIV)
#endif

//to use most bound center
//~ #define SUBCOM(i) Pdat.Pos[SubCat.PSubArr[i][0]]  
//to use core com (or minpot in case of subfind)
#define SUBCOM(i) SubCat.Property[i].CoM

//~ #define SUBFIND_DIR "/home/kambrain/data/6702DM/subcatS"   //define this properly to apply on subfind, 
															//undef for HBT

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(HBTInt Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif


HBTInt *PartHost;
SUBCATALOGUE SubCat;
CATALOGUE Cat;

typedef struct 
{
	HBTInt nbin;
	HBTInt flag_badbins;//1 means badbin,-1 means bad host,0 normal,-2 means rcen<rmin
	HBTReal rmax;
	HBTReal rmin;//min radius of radial log-bins; equals to max(Factor_RMIN*rmax,SofteningHalo); note this is not saved.
	HBTReal rvir;
	HBTReal req_bk_1;//set to 0 when no bin,to rmax when no back
	HBTReal req_all_1;//set to 0 when no bin,to rmax when no back+other
	HBTReal req_bk_02;
	HBTReal req_all_02;
	HBTReal rtidal;//set to rvir when cannot determine (e.g.,mainsubs/quasi-subs and subhalos inside a bad-host )
	HBTReal rcen;
	
	HBTInt *n_this;
	HBTInt *n_other;
	HBTInt *n_back;
	HBTReal *r_this;
	HBTReal *r_other;
	HBTReal *r_back;
} SUBPROFILE;

void sub_radius_eq_local(SUBPROFILE *subprof);//consider improving accuracy by interpolation or fitting
void sub_radius_eq_global(SUBPROFILE *subprof);//consider improving accuracy by interpolation or fitting
void sub_radius_eq_aqua(SUBPROFILE *subprof,HBTInt SubLen);
void sub_bin(HBTInt subid,SUBPROFILE *subprof);
HBTInt fill_subprof(SUBPROFILE *subprof,HBTReal *rp_this,HBTReal *rp_other,HBTReal *rp_back,
									     HBTInt np_this,HBTInt np_other,HBTInt np_back);
HBTInt histr(HBTReal *rdata,HBTInt ndata,HBTReal rmax,HBTInt *count,HBTReal *rbin,HBTInt nbin);
HBTInt histr_log(HBTReal *rdata,HBTInt ndata,HBTReal rmin,HBTReal rmax,HBTInt *count,HBTReal *rbin,HBTInt nbin);
void allocate_subprof(SUBPROFILE *subprof,HBTInt nbin);
void free_subprof(SUBPROFILE *subprof);
void save_subs_size(SUBPROFILE *subprof,HBTInt Nsubs,HBTInt Nsnap,char *outputdir);
void save_subs_prof(SUBPROFILE *subprof,HBTInt Nsubs,HBTInt Nsnap,char *outputdir);

int main(int argc, char **argv)
{
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	
HBTInt Nsnap=99,i;
HBTInt subid,grpid,hostsubid,hostnbin,snbin,*sn_host,*n_host,icen;
HBTReal dmr,vratio,hostrmax,hostrmin;
double hostspan;
SUBPROFILE *subprof;

if(argc!=2)
{
  printf("Usage: %s [Nsnap]    ; now exit\n", argv[0]);
  exit(1);
}
Nsnap=atoi(argv[1]);

logfile=stdout;
sprintf(outputdir,"%s/profile/logbin",SUBCAT_DIR);
mkdir(outputdir,0755);
#ifdef AQUA_SIZE
sprintf(outputdir,"%s/profile/logbin/aqua",SUBCAT_DIR);		
#endif
mkdir(outputdir,0755);

load_sub_catalogue(Nsnap,&SubCat,inputdir);
load_particle_data(Nsnap,SNAPSHOT_DIR);
fill_PIDHash();
fresh_ID2Index(&SubCat,-2);
free_PIDHash();

PartHost=mymalloc(sizeof(HBTInt)*NP_DM);
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
	#ifdef AQUA_SIZE
	sub_radius_eq_aqua(subprof+subid,SubCat.SubLen[subid]);//calculating req as done in Aquarius
	#else
	sub_radius_eq_global(subprof+subid);//calculating req
	#endif
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
			hostrmin=subprof[hostsubid].rmin;
			hostspan=log(hostrmax/hostrmin);
			vratio=(1-exp(-3.*hostspan/hostnbin));
			//~ n_host=subprof[hostsubid].n_this;
			n_host=mymalloc(sizeof(HBTInt)*hostnbin);//all as tidal
			sn_host=mymalloc(sizeof(HBTInt)*hostnbin);
			snbin=0;
			for(i=0;i<hostnbin;i++)
			{
				n_host[i]=subprof[hostsubid].n_this[i]+subprof[hostsubid].n_other[i]+subprof[hostsubid].n_back[i];//all as tidal
				snbin+=n_host[i];
				sn_host[i]=snbin;
			}
			for(subid=hostsubid+1;subid<hostsubid+SubCat.GrpLen_Sub[grpid];subid++)
			{
				//~ icen=floor(subprof[subid].rcen/hostrmax*hostnbin);
				icen=floor(log(subprof[subid].rcen/hostrmin)/hostspan*hostnbin);
				if(icen>=hostnbin)//a far-away sub
				{
				subprof[subid].flag_badbins=2;//more than outer	
				icen=hostnbin-1;
				dmr=0.;
				}
				else if(icen<0)
				{	
				subprof[subid].flag_badbins=-2;//more than inner
				icen=0;
				dmr=3.;
				}
				else
				{
				//~ vratio=(HBTReal)icen/(icen+1);
				//~ vratio=vratio*vratio*vratio;
				dmr=3.*n_host[icen]/vratio/sn_host[icen]*(1-exp(-3.*hostspan*(icen+1)/hostnbin));
				}
				if(dmr>=TIDAL_MAX)//probably a badly sampled or a quite loose host (badbin)
				{
				subprof[subid].flag_badbins=1;
				subprof[subid].rtidal=subprof[subid].rvir;
				}
				else
				subprof[subid].rtidal=subprof[subid].rcen*pow(SubCat.SubLen[subid]/(TIDAL_MAX-dmr)/sn_host[icen],1.0/3);
			}
			free(n_host);//all as tidal
			free(sn_host);
		}
		else//bad_host
		{
		  for(subid=hostsubid+1;subid<hostsubid+SubCat.GrpLen_Sub[grpid];subid++)//2013-05-07 22:14:48 bug fix: added this line.
		  {
			subprof[subid].flag_badbins=-1;
			subprof[subid].rtidal=subprof[subid].rvir;
		  }
		}
	}
}
for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++)//quasi
	subprof[subid].rtidal=subprof[subid].rvir;

load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);	
save_subs_size(subprof,SubCat.Nsubs,Nsnap,outputdir);
save_subs_prof(subprof,SubCat.Nsubs,Nsnap,outputdir);
free_catalogue(&Cat);
for(subid=0;subid<SubCat.Nsubs;subid++)
	free_subprof(subprof+subid);
myfree(subprof);

myfree(PartHost);	
erase_sub_catalogue(&SubCat);
free_particle_data();

return 0;
}

void sub_radius_eq_aqua(SUBPROFILE *subprof,HBTInt SubLen)//consider improving accuracy by interpolation or fitting
{/*the radius where subdens equals to global  background
* (background averaged inside rmax,not inside rtidal,
* since inside rtidal the background level can be too low)*
* * calculate background inside a radius r with Mtotal(<r)=Msub*/
	HBTReal *diff_dens_back,*diff_dens_all,av_dens_back,av_dens_all,*vol;
	HBTReal rmin,rmax,vratbase;
	double span;
	HBTInt i,nbin,nbin_back,sum_this,sum_back,sum_other,binmax_sub;
	if(nbin=subprof->nbin)
	{
		sum_this=0;sum_back=0;sum_other=0;
		for(i=0;i<nbin;i++)
		{
			if(sum_this+sum_back+sum_other>=SubLen)
			break;
			sum_this+=subprof->n_this[i];
			sum_back+=subprof->n_back[i];
			sum_other+=subprof->n_other[i];
		}
		nbin_back=i;
		if(sum_back+sum_other)
		{
			rmin=subprof->rmin;
			rmax=subprof->rmax;
			span=log(rmax/rmin);
			vratbase=(exp(3.*span/nbin)-1)/(exp(3.*nbin_back*(HBTReal)span/nbin)-1);
			vol=mymalloc(sizeof(HBTReal)*nbin);
			for(i=0;i<nbin;i++)
				vol[i]=3.*i*span/nbin;
				
			diff_dens_all=mymalloc(sizeof(HBTReal)*nbin);	
			av_dens_all=(HBTReal)(sum_back+sum_other)*vratbase;
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=logf((HBTReal)subprof->n_this[i]/av_dens_all)-vol[i]; //  logf(/exp(vol[i]))=-vol[i]
			binmax_sub=max_of_vec(diff_dens_all,nbin);
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(diff_dens_all[i]);
			subprof->req_all_1=subprof->r_this[binmax_sub+min_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			av_dens_all*=0.02;
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(logf((HBTReal)subprof->n_this[i]/av_dens_all)-vol[i]);
			subprof->req_all_02=subprof->r_this[binmax_sub+min_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			free(diff_dens_all);
			
			if(sum_back)
			{
				diff_dens_back=mymalloc(sizeof(HBTReal)*nbin);
				av_dens_back=(HBTReal)sum_back*vratbase;
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((HBTReal)subprof->n_this[i]/av_dens_back)-vol[i]);
				subprof->req_bk_1=subprof->r_this[binmax_sub+min_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
				av_dens_back*=0.02;
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((HBTReal)subprof->n_this[i]/av_dens_back)-vol[i]);
				subprof->req_bk_02=subprof->r_this[binmax_sub+min_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
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
void sub_radius_eq_global(SUBPROFILE *subprof)//consider improving accuracy by interpolation or fitting
{/*the radius where subdens equals to global  background
* (background averaged inside rmax,not inside rtidal,
* since inside rtidal the background level can be too low)*/
	HBTReal *diff_dens_back,*diff_dens_all,av_dens_back,av_dens_all,*vol;
	HBTReal rmin,rmax,vratbase;
	double span;
	HBTInt i,nbin,sum_back,sum_other,binmax_sub;
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
			rmin=subprof->rmin;
			rmax=subprof->rmax;
			span=log(rmax/rmin);
			vratbase=(exp(3.*span/nbin)-1)/(exp(3.*span)-1);
			vol=mymalloc(sizeof(HBTReal)*nbin);
			for(i=0;i<nbin;i++)
				vol[i]=3.*i*span/nbin;
				
			diff_dens_all=mymalloc(sizeof(HBTReal)*nbin);	
			av_dens_all=(HBTReal)(sum_back+sum_other)*vratbase;
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=logf((HBTReal)subprof->n_this[i]/av_dens_all)-vol[i]; //  logf(/exp(vol[i]))=-vol[i]
			binmax_sub=max_of_vec(diff_dens_all,nbin);
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(diff_dens_all[i]);
			subprof->req_all_1=subprof->r_this[binmax_sub+min_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			av_dens_all*=0.02;
			for(i=0;i<nbin;i++)
				diff_dens_all[i]=fabs(logf((HBTReal)subprof->n_this[i]/av_dens_all)-vol[i]);
			subprof->req_all_02=subprof->r_this[binmax_sub+min_of_vec(diff_dens_all+binmax_sub,nbin-binmax_sub)];
			free(diff_dens_all);
			
			if(sum_back)
			{
				diff_dens_back=mymalloc(sizeof(HBTReal)*nbin);
				av_dens_back=(HBTReal)sum_back*vratbase;
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((HBTReal)subprof->n_this[i]/av_dens_back)-vol[i]);
				subprof->req_bk_1=subprof->r_this[binmax_sub+min_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
				av_dens_back*=0.02;
				for(i=0;i<nbin;i++)
					diff_dens_back[i]=fabs(logf((HBTReal)subprof->n_this[i]/av_dens_back)-vol[i]);
				subprof->req_bk_02=subprof->r_this[binmax_sub+min_of_vec(diff_dens_back+binmax_sub,nbin-binmax_sub)];
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
	HBTReal *diff_dens_back,*diff_dens_all;
	HBTInt i;
	if(subprof->nbin)
	{
	diff_dens_back=mymalloc(sizeof(HBTReal)*subprof->nbin);
	diff_dens_all=mymalloc(sizeof(HBTReal)*subprof->nbin);

	for(i=0;i<subprof->nbin;i++)
	{
		diff_dens_back[i]=fabs(subprof->n_this[i]-subprof->n_back[i])/(HBTReal)subprof->n_this[i];
		diff_dens_all[i]=fabs(subprof->n_this[i]-subprof->n_back[i]-subprof->n_other[i])/(HBTReal)subprof->n_this[i];
	}
	subprof->req_bk_1=subprof->r_this[min_of_vec(diff_dens_back,subprof->nbin)];
	subprof->req_all_1=subprof->r_this[min_of_vec(diff_dens_all,subprof->nbin)];
	
	for(i=0;i<subprof->nbin;i++)
	{
		diff_dens_back[i]=fabs(subprof->n_this[i]-0.02*subprof->n_back[i])/(HBTReal)subprof->n_this[i];
		diff_dens_all[i]=fabs(subprof->n_this[i]-0.02*(subprof->n_back[i]+subprof->n_other[i]))/(HBTReal)subprof->n_this[i];
	}
	subprof->req_bk_02=subprof->r_this[min_of_vec(diff_dens_back,subprof->nbin)];
	subprof->req_all_02=subprof->r_this[min_of_vec(diff_dens_all,subprof->nbin)];
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

void sub_bin(HBTInt subid,SUBPROFILE *subprof)
{
	HBTReal *cen,rmax,dr;
	HBTReal *rp_this,*rp_other,*rp_back;
	HBTInt    np_this,np_other,np_back,
		   np_this_max,np_other_max,np_back_max;
	HBTInt i,j,k,pid;
	HBTInt subbox_grid[3][2];
	
	if(SubCat.SubLen[subid])
	{
	rp_this=mymalloc(sizeof(HBTReal)*SubCat.SubLen[subid]);
	rp_other=mymalloc(sizeof(HBTReal)*SubCat.SubLen[subid]);
	rp_back=mymalloc(sizeof(HBTReal)*SubCat.SubLen[subid]);
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
							rp_this=realloc(rp_this,sizeof(HBTReal)*np_this_max);
							}
							rp_this[np_this]=dr;
							np_this++;
						}
						else if(PartHost[pid]==-1)
						{
							if(np_back>=np_back_max)
							{
							np_back_max*=2;	
							rp_back=realloc(rp_back,sizeof(HBTReal)*np_back_max);
							}
							rp_back[np_back]=dr;
							np_back++;
						}
						else
						{
							if(np_other>=np_other_max)
							{
							np_other_max*=2;	
							rp_other=realloc(rp_other,sizeof(HBTReal)*np_other_max);
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
	subprof->rmin=0;
	}
}
HBTInt fill_subprof(SUBPROFILE *subprof,HBTReal *rp_this,HBTReal *rp_other,HBTReal *rp_back,
									     HBTInt np_this,HBTInt np_other,HBTInt np_back)
{
	HBTReal rmax,rmin;
	HBTInt i,nbin,flag_rebin;

	nbin=NBIN_MAX;
	subprof->n_this=mymalloc(sizeof(HBTInt)*nbin);
	subprof->r_this=mymalloc(sizeof(HBTReal)*nbin);
	rmax=rp_this[max_of_vec(rp_this,np_this)];
	subprof->rmax=rmax;
	rmin=rmax*Factor_RMIN;
	if(rmin<SofteningHalo)
		rmin=SofteningHalo;
	subprof->rmin=rmin;
	flag_rebin=1;nbin++;
	while(flag_rebin&&nbin>NBIN_MIN)
	{
	nbin--;
	flag_rebin=histr_log(rp_this,np_this,rmin,rmax,subprof->n_this,subprof->r_this,nbin);
	}
	if(NBIN_MIN==nbin)
	{
	free(subprof->n_this);subprof->n_this=NULL;
	free(subprof->r_this);subprof->r_this=NULL;
	subprof->nbin=0;
	return 0;
	}
	subprof->nbin=nbin;
	subprof->n_this=realloc(subprof->n_this,sizeof(HBTInt)*nbin);
	subprof->r_this=realloc(subprof->r_this,sizeof(HBTReal)*nbin);
	subprof->n_other=mymalloc(sizeof(HBTInt)*nbin);
	subprof->r_other=mymalloc(sizeof(HBTReal)*nbin);
	subprof->n_back=mymalloc(sizeof(HBTInt)*nbin);
	subprof->r_back=mymalloc(sizeof(HBTReal)*nbin);
	histr_log(rp_other,np_other,rmin,rmax,subprof->n_other,subprof->r_other,nbin);
	histr_log(rp_back,np_back,rmin,rmax,subprof->n_back,subprof->r_back,nbin);
	return nbin;
}

HBTInt histr(HBTReal *rdata,HBTInt ndata,HBTReal rmax,HBTInt *count,HBTReal *rbin,HBTInt nbin)
{
	//return 1 if have empty bins,0 otherwise
	HBTInt i,dr_bin,flag_emptybin;
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
HBTInt histr_log(HBTReal *rdata,HBTInt ndata,HBTReal rmin,HBTReal rmax,HBTInt *count,HBTReal *rbin,HBTInt nbin)
{
	//return 1 if have empty bins,0 otherwise
	HBTInt i,dr_bin,flag_emptybin;
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
void allocate_subprof(SUBPROFILE *subprof,HBTInt nbin)
{
	subprof->nbin=nbin;
	if(nbin)
	{
		subprof->n_this=mymalloc(sizeof(HBTInt)*nbin);
		subprof->n_other=mymalloc(sizeof(HBTInt)*nbin);
		subprof->n_back=mymalloc(sizeof(HBTInt)*nbin);
		subprof->r_this=mymalloc(sizeof(HBTReal)*nbin);
		subprof->r_other=mymalloc(sizeof(HBTReal)*nbin);
		subprof->r_back=mymalloc(sizeof(HBTReal)*nbin);
	}
}

void free_subprof(SUBPROFILE *subprof)
{
	if(subprof->nbin)
	{
		myfree(subprof->n_this);
		myfree(subprof->n_other);
		myfree(subprof->n_back);
		myfree(subprof->r_this);
		myfree(subprof->r_other);
		myfree(subprof->r_back);
	}
}

//~ void save_subs_size(SUBPROFILE *subprof,HBTInt Nsubs,HBTInt Nsnap,char *outputdir)
//~ {
	//~ HBTInt i;
	//~ char buf[1024];
	//~ FILE *fp;
	//~ sprintf(buf,"%s/sub_size_%03d",outputdir,Nsnap);
	//~ myfopen(fp,buf,"w");
	//~ for(i=0;i<Nsubs;i++)
	//~ {
		//~ fwrite(&subprof[i].nbin,sizeof(HBTInt),1,fp);
		//~ fwrite(&subprof[i].flag_badbins,sizeof(HBTInt),1,fp);
		//~ fwrite(&subprof[i].rmax,sizeof(HBTReal),1,fp);
		//~ fwrite(&subprof[i].rvir,sizeof(HBTReal),1,fp);
		//~ fwrite(&subprof[i].req_bk_1,sizeof(HBTReal),1,fp);
		//~ fwrite(&subprof[i].req_all_1,sizeof(HBTReal),1,fp);
		//~ fwrite(&subprof[i].req_bk_02,sizeof(HBTReal),1,fp);
		//~ fwrite(&subprof[i].req_all_02,sizeof(HBTReal),1,fp);
		//~ fwrite(&subprof[i].rtidal,sizeof(HBTReal),1,fp);
		//~ fwrite(&subprof[i].rcen,sizeof(HBTReal),1,fp);
	//~ }
	//~ fclose(fp);
//~ }
//~ void save_subs_prof(SUBPROFILE *subprof,HBTInt Nsubs,HBTInt Nsnap,char *outputdir)
//~ {
	//~ HBTInt i,nbin;
	//~ char buf[1024];
	//~ FILE *fp;
	//~ sprintf(buf,"%s/sub_prof_%03d",outputdir,Nsnap);
	//~ myfopen(fp,buf,"w");
	//~ for(i=0;i<Nsubs;i++)
	//~ {
		//~ nbin=subprof[i].nbin;
		//~ fwrite(subprof[i].n_this,sizeof(HBTInt),nbin,fp);
		//~ fwrite(subprof[i].n_other,sizeof(HBTInt),nbin,fp);
		//~ fwrite(subprof[i].n_back,sizeof(HBTInt),nbin,fp);
		//~ fwrite(subprof[i].r_this,sizeof(HBTReal),nbin,fp);
		//~ fwrite(subprof[i].r_other,sizeof(HBTReal),nbin,fp);
		//~ fwrite(subprof[i].r_back,sizeof(HBTReal),nbin,fp);
	//~ }
	//~ fclose(fp);
//~ }
void save_subs_size(SUBPROFILE *subprof,HBTInt Nsubs,HBTInt Nsnap,char *outputdir)
{
	HBTInt grpid,i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/sat_size_%03d.%d",outputdir,(int)Nsnap,TIDAL_MAX);
	myfopen(fp,buf,"w");
	for(grpid=0;grpid<GRPID_MAX;grpid++)
	{
		if(Cat.Len[grpid]<GRPLEN_MIN) break;
		for(i=SubCat.GrpOffset_Sub[grpid]+1;i<SubCat.GrpOffset_Sub[grpid]+SubCat.GrpLen_Sub[grpid];i++)
		{
			if(SubCat.SubLen[i]<SUBLEN_MIN) break;
			fwrite(&i,sizeof(HBTInt),1,fp);
			fwrite(&subprof[i].nbin,sizeof(HBTInt),1,fp);
			fwrite(&subprof[i].flag_badbins,sizeof(HBTInt),1,fp);
			fwrite(&subprof[i].rmax,sizeof(HBTReal),1,fp);
			fwrite(&subprof[i].rvir,sizeof(HBTReal),1,fp);
			fwrite(&subprof[i].req_bk_1,sizeof(HBTReal),1,fp);
			fwrite(&subprof[i].req_all_1,sizeof(HBTReal),1,fp);
			fwrite(&subprof[i].req_bk_02,sizeof(HBTReal),1,fp);
			fwrite(&subprof[i].req_all_02,sizeof(HBTReal),1,fp);
			fwrite(&subprof[i].rtidal,sizeof(HBTReal),1,fp);
			fwrite(&subprof[i].rcen,sizeof(HBTReal),1,fp);
		}
	}
	fclose(fp);
	sprintf(buf,"%s/main_size_%03d.%d",outputdir,(int)Nsnap,TIDAL_MAX);
	myfopen(fp,buf,"w");
	for(grpid=0;grpid<GRPID_MAX;grpid++)
	{
		if(Cat.Len[grpid]<GRPLEN_MIN) break;
		if(0==SubCat.GrpLen_Sub[grpid]) continue;  //this for subfind catalogue where offset do not increase when halo has no sub
		i=SubCat.GrpOffset_Sub[grpid];
		if(SubCat.SubLen[i]<SUBLEN_MIN) break;
		fwrite(&i,sizeof(HBTInt),1,fp);
		fwrite(&subprof[i].nbin,sizeof(HBTInt),1,fp);
		fwrite(&subprof[i].flag_badbins,sizeof(HBTInt),1,fp);
		fwrite(&subprof[i].rmax,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].rvir,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_bk_1,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_all_1,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_bk_02,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_all_02,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].rtidal,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].rcen,sizeof(HBTReal),1,fp);
	}
	for(i=SubCat.Nsubs-SubCat.NQuasi;i<SubCat.Nsubs;i++)//quasi
	{
		if(SubCat.SubLen[i]<SUBLEN_MIN) break;
		fwrite(&i,sizeof(HBTInt),1,fp);
		fwrite(&subprof[i].nbin,sizeof(HBTInt),1,fp);
		fwrite(&subprof[i].flag_badbins,sizeof(HBTInt),1,fp);
		fwrite(&subprof[i].rmax,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].rvir,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_bk_1,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_all_1,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_bk_02,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].req_all_02,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].rtidal,sizeof(HBTReal),1,fp);
		fwrite(&subprof[i].rcen,sizeof(HBTReal),1,fp);
	}
	fclose(fp);
}
void save_subs_prof(SUBPROFILE *subprof,HBTInt Nsubs,HBTInt Nsnap,char *outputdir)
{
	HBTInt grpid,i,nbin;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/sat_prof_%03d.%d",outputdir,(int)Nsnap,TIDAL_MAX);
	myfopen(fp,buf,"w");
	for(grpid=0;grpid<GRPID_MAX;grpid++)
	{
		if(Cat.Len[grpid]<GRPLEN_MIN) break;
		for(i=SubCat.GrpOffset_Sub[grpid]+1;i<SubCat.GrpOffset_Sub[grpid]+SubCat.GrpLen_Sub[grpid];i++)
		{
		if(SubCat.SubLen[i]<SUBLEN_MIN) break;
		nbin=subprof[i].nbin;
		fwrite(subprof[i].n_this,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].n_other,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].n_back,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].r_this,sizeof(HBTReal),nbin,fp);
		fwrite(subprof[i].r_other,sizeof(HBTReal),nbin,fp);
		fwrite(subprof[i].r_back,sizeof(HBTReal),nbin,fp);
		}
	}
	fclose(fp);
	sprintf(buf,"%s/main_prof_%03d.%d",outputdir,(int)Nsnap,TIDAL_MAX);
	myfopen(fp,buf,"w");
	for(grpid=0;grpid<GRPID_MAX;grpid++)
	{
		if(Cat.Len[grpid]<GRPLEN_MIN) break;
		if(0==SubCat.GrpLen_Sub[grpid]) continue;  //this for subfind catalogue where offset do not increase when halo has no sub
		i=SubCat.GrpOffset_Sub[grpid];
		if(SubCat.SubLen[i]<SUBLEN_MIN) break;
		nbin=subprof[i].nbin;
		fwrite(subprof[i].n_this,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].n_other,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].n_back,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].r_this,sizeof(HBTReal),nbin,fp);
		fwrite(subprof[i].r_other,sizeof(HBTReal),nbin,fp);
		fwrite(subprof[i].r_back,sizeof(HBTReal),nbin,fp);
	}
	for(i=SubCat.Nsubs-SubCat.NQuasi;i<SubCat.Nsubs;i++)//quasi
	{
		if(SubCat.SubLen[i]<SUBLEN_MIN) break;
		nbin=subprof[i].nbin;
		fwrite(subprof[i].n_this,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].n_other,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].n_back,sizeof(HBTInt),nbin,fp);
		fwrite(subprof[i].r_this,sizeof(HBTReal),nbin,fp);
		fwrite(subprof[i].r_other,sizeof(HBTReal),nbin,fp);
		fwrite(subprof[i].r_back,sizeof(HBTReal),nbin,fp);
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

/* changelog:
 * 2013-05-07 20:43:12 : adapted to the new datatypes and linklist functions
 * 2013-05-07 22:26:25 : bug fix around line 213.
 * 2013-05-08 01:12:22 : bug fix: previously forgot to output quasi's.
 * */

