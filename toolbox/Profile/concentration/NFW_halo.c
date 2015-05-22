#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "datatypes.h"

#define PRINT_PROF

#define NUM_PARAM 2

#define TOL_FIT 1e-6
#define MAX_ITER 500

#define FIT_FOF  //fit fof halo or spherical halo
//#define FIT_VIR  //fit only density inside virial radius

//#define PSEUDO_LOGFIT  //using sigma=rho as weight to mimic fiting in logspace
#define LOGFIT //fit in logspace, no effect if PSEUDO_LOGFIT is defined
#define POISSON_ERR //use poisson error as weight, no effect if PSEUDO_LOGFIT is defined

//undef this only if you have resonable error for your data
//#define ERR_EST_FIT_ALWAYS  
							/*always give error estimation for parameters from fitting,
							 *define this if your data's error is not serious 
							 * but just meant to be a weight, or if you do not 
							 * have error from data.
							 * otherwise for good fits (small chisq/dof) 
							 * error is derived only from error propagation
							 * (because in this case the provided data's error 
							 * outweight the uncertainty from fitting)*/
							
#ifdef PSEUDO_LOGFIT 
#undef LOGFIT
#undef POISSON_ERR
#endif

#ifdef LOGFIT
#include "lnNFW_fun.c"
#else
#include "NFW_fun.c"
#endif

#define Factor_RMIN 1e-2
#define mymalloc malloc
#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ fprintf(stderr,"Error opening file '%s'\n",filename);	fflush(stderr); exit(1);	}

#define NFW_SUCCESS 0
#define NFW_BADBIN 1
#define NFW_FAIL 2

typedef struct 
{
	int nbin;
	float rmax;
	float rmin;
	float rvir;//comoving tophat rvir estimated from fof mass
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

void load_halos_prof(HALOPROFILE *haloprof,int Ngrps,int Nsnap);
void free_haloprof(HALOPROFILE *haloprof);
HBTInt read_ngroups(int Nsnap);
int fit_halo_prof(HALOPROFILE *haloprof,double par[2],double err[2],double err_M[3]);

int main (int argc,char **argv)
{
double par[2],err[2],err_M[3];
int status;
HALOPROFILE *haloprof;	
int Nsnap,Ngroups,grpid;

Nsnap=atoi(argv[1]);
grpid=atoi(argv[2]);

Ngroups=read_ngroups(Nsnap);
haloprof=mymalloc(sizeof(HALOPROFILE)*Ngroups);
load_halos_prof(haloprof,Ngroups,Nsnap);

status=fit_halo_prof(haloprof+grpid,par,err,err_M);

for(grpid=0;grpid<Ngroups;grpid++)
	free_haloprof(haloprof+grpid);
free(haloprof);

return 0;
}

void load_halos_prof(HALOPROFILE *haloprof,int Ngrps,int Nsnap)
{
	int i,nbin;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngrps;i++)
	{
		fread(&haloprof[i].nbin,sizeof(int),1,fp);
		fread(&haloprof[i].rmax,sizeof(float),1,fp);
		fread(&haloprof[i].rvir,sizeof(float),1,fp);
		haloprof[i].rmin=Factor_RMIN*haloprof[i].rmax;
		if(haloprof[i].rmin<SofteningHalo) haloprof[i].rmin=SofteningHalo;
		fseek(fp,10*4L,SEEK_CUR);
		fread(&haloprof[i].mass,sizeof(int),1,fp);
		fread(haloprof[i].Mvir,sizeof(int),3,fp);
		fread(haloprof[i].Rvir,sizeof(float),3,fp);
		fread(haloprof[i].flag_badvir,sizeof(int),3,fp);
		fread(&haloprof[i].flag_fakehalo,sizeof(int),1,fp);
	}
	fclose(fp);
	sprintf(buf,"%s/profile/logbin/halo_prof_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngrps;i++)
	{
		nbin=haloprof[i].nbin;
		haloprof[i].n_this=mymalloc(sizeof(int)*nbin);
		haloprof[i].n_other=mymalloc(sizeof(int)*nbin);
		haloprof[i].n_back=mymalloc(sizeof(int)*nbin);
		haloprof[i].r_this=mymalloc(sizeof(float)*nbin);
		haloprof[i].r_other=mymalloc(sizeof(float)*nbin);
		haloprof[i].r_back=mymalloc(sizeof(float)*nbin);		
		fread(haloprof[i].n_this,sizeof(int),nbin,fp);
		fread(haloprof[i].n_other,sizeof(int),nbin,fp);
		fread(haloprof[i].n_back,sizeof(int),nbin,fp);
		fread(haloprof[i].r_this,sizeof(float),nbin,fp);
		fread(haloprof[i].r_other,sizeof(float),nbin,fp);
		fread(haloprof[i].r_back,sizeof(float),nbin,fp);
	}
	fclose(fp);
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
HBTInt read_ngroups(int Nsnap)
{
	HBTInt ngroups;
	char buf[1024];
	FILE *fp;
	
	 sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, Nsnap);
 	 myfopen(fp,buf,"r");
   fread(&ngroups, sizeof(HBTInt), 1, fp);
   fclose(fp);
   return ngroups;
}

double * volume_logbin(float rmin,float rmax,int nbin)
{
	int i;
	double *vol,vl,vr,vbin;
	const double c=4./3*M_PI;
	vol=mymalloc(sizeof(double)*nbin);
	vbin=exp((log(rmax)-log(rmin))/nbin*3.);
	vl=c*pow(rmin,3.);
	for(i=0;i<nbin;i++)
	{
		vr=vl*vbin;
		vol[i]=vr-vl;
		vl=vr;
	}
	return vol;
}

int set_solver_data(HALOPROFILE *haloprof,gsl_multifit_function_fdf *f)
{
double *r, *y, *sigma, *vol;
int na,nbin,i;
struct data *d;

nbin=haloprof->nbin;
if(nbin<=NUM_PARAM) 
{
	printf("Ndata=%d\n",nbin);
	return 0;
}

vol=volume_logbin(haloprof->rmin,haloprof->rmax,nbin);

#ifdef FIT_VIR
nbin=floor(nbin*log(haloprof->Rvir[0]/haloprof->rmin)
				/log(haloprof->rmax/haloprof->rmin))+1;
if(nbin>haloprof->nbin) nbin=haloprof->nbin;				
#endif
printf("Ndata=%d\n",nbin);
if(nbin<=NUM_PARAM) return 0;

/* This is the data to be fitted */
d=mymalloc(sizeof(struct data));
r=mymalloc(sizeof(double)*nbin);
y=mymalloc(sizeof(double)*nbin);
sigma=mymalloc(sizeof(double)*nbin);

for (i = 0; i < nbin; i++)
 {
   #ifdef FIT_FOF
   na=haloprof->n_this[i];
   r[i]=haloprof->r_this[i];
   #else
   na=haloprof->n_this[i]+haloprof->n_other[i]+haloprof->n_back[i];
   r[i]=haloprof->n_this[i]*haloprof->r_this[i]+
   	  haloprof->n_other[i]*haloprof->r_other[i]+
	  haloprof->n_back[i]*haloprof->r_back[i];
   r[i]=r[i]/na;
   #endif
   y[i]=na/vol[i];
   #ifdef PSEUDO_LOGFIT
	   sigma[i]=y[i];
   #else
  	   #ifdef POISSON_ERR
		   sigma[i]=sqrt(na)/vol[i];
	   #else
		   sigma[i]=1;
	   #endif
	   #ifdef LOGFIT
		   #ifdef POISSON_ERR
		   sigma[i]=1./sqrt(na);// dlog(y)=dy/y
		   #endif
		   y[i]=log(y[i]);
	   #endif  
   #endif
   //~ printf ("data: %u %g %g %g\n", i,r[i], y[i], sigma[i]);
 };
 free(vol);
 
d->nbin=nbin;
d->r=r;
d->rho=y;
d->sigma=sigma; 

#ifdef LOGFIT
f->f = &lnNFW_f;
f->df = &lnNFW_df;
f->fdf = &lnNFW_fdf;
#else
f->f = &NFW_f;
f->df = &NFW_df;
f->fdf = &NFW_fdf;
#endif
f->n = nbin;
f->p = NUM_PARAM;
f->params = d;

return 1;
}

void free_solver_data(gsl_multifit_function_fdf *f)
{
struct data *d;
d=f->params;
free(d->r);
free(d->rho);
free(d->sigma);
free(d);
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
printf ("iter: %ld x = % 15.8f % 15.8f "
	   "|f(x)| = %g\n",	   iter,
	   gsl_vector_get (s->x, 0), 
	   gsl_vector_get (s->x, 1),
	   gsl_blas_dnrm2 (s->f));
}

int fit_halo_prof(HALOPROFILE *haloprof,double par[2],double err[2],double err_M[3])
{	
gsl_multifit_function_fdf f;
const gsl_multifit_fdfsolver_type *T;
gsl_multifit_fdfsolver *s;
gsl_matrix *covar;
int status;
unsigned int i,iter = 0;
double x_init[2];// = { 0., 100.0};
//initial guess, using c=5;
x_init[0]=0.;
x_init[1]=haloprof->rvir/5.;
gsl_vector_view x = gsl_vector_view_array (x_init, NUM_PARAM);

status=set_solver_data(haloprof,&f);
if(status==0)
{
 par[0]=par[1]=0.;
 err[0]=err[1]=0.;
 err_M[0]=err_M[1]=err_M[2]=0.;	
 return NFW_BADBIN;
}

#ifdef PRINT_PROF
struct data *d=f.params;
printf("-------------------------------------------\n")
printf("#r, rho, sigma\n");
int i;
for(i=0;i<d->nbin;i++)
 printf("%g, %g, %g\n", d->r[i], d->rho[i], d->sigma[i]);
printf("-------------------------------------------\n")
#endif

covar = gsl_matrix_alloc (NUM_PARAM, NUM_PARAM);
T = gsl_multifit_fdfsolver_lmsder;
s = gsl_multifit_fdfsolver_alloc (T, f.n, f.p);
gsl_multifit_fdfsolver_set (s, &f, &x.vector);
print_state (iter, s);
do
 {
   iter++;
   status = gsl_multifit_fdfsolver_iterate (s);
   printf ("status = %s\n", gsl_strerror (status));
   print_state (iter, s);
   if (status)
	 break;
   status = gsl_multifit_test_delta (s->dx, s->x, TOL_FIT, TOL_FIT);
 }
while (status == GSL_CONTINUE && iter < MAX_ITER);
gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
{ 
 double chi = gsl_blas_dnrm2(s->f);//2-norm of the residual function
 double dof = f.n-f.p;
 #ifdef ERR_EST_FIT_ALWAYS
 double c = chi / sqrt(dof); 
 #else
 double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
 #endif
 printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
#ifdef LOGFIT
 par[0]=exp(FIT(0));
 err[0]=par[0]*c*ERR(0);
#else
 par[0]=FIT(0);err[0]=c*ERR(0);
#endif
 par[1]=FIT(1);err[1]=c*ERR(1);
 printf ("rhos = %.3e +/- %.2e\n", par[0]*MP_DM, err[0]*MP_DM);
 printf ("rs   = %.3e +/- %2e\n", par[1], err[1]);
}
#define CONCEN(i) (haloprof->Rvir[i]/par[1])
#define MFIT(i) (4*M_PI*par[1]*par[1]*par[1]*par[0]*(log(1+CONCEN(i))-CONCEN(i)/(1+CONCEN(i))))
for(i=0;i<3;i++)
{
	if(haloprof->flag_badvir[i])
	err_M[i]=fabs(MFIT(i)/(float)(haloprof->mass)-1.);
	else
	err_M[i]=fabs(MFIT(i)/(float)(haloprof->Mvir[i])-1.);
}
printf("Precision of Fitted Virial Mass:%.2g, %.2g, %.2g\n",err_M[0],err_M[1],err_M[2]);
printf("concentration= %g, %g, %g\n",CONCEN(0),CONCEN(1),CONCEN(2));
printf("Mvir(TH,C200,B200): %g, %g, %g\n",(float)haloprof->Mvir[0]*MP_DM,(float)haloprof->Mvir[1]*MP_DM,(float)haloprof->Mvir[2]*MP_DM);
printf ("iter=%u, status = %s\n", iter, gsl_strerror (status));

free_solver_data(&f);
gsl_multifit_fdfsolver_free (s);			
gsl_matrix_free (covar);
if(status==GSL_SUCCESS)
return NFW_SUCCESS;
else
return NFW_FAIL;
}
