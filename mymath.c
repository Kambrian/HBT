#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

HBTInt Dmax_of_vec(double *vec,HBTInt Len)//what if the max is not unique?,return first max
{HBTInt tmp,i;
	tmp=0;
	for(i=0;i<Len;i++)
	{
		if(vec[tmp]<vec[i])tmp=i;
	}
	return tmp;//max ID
}
HBTInt Fmax_of_vec(float *vec,HBTInt Len)//what if the max is not unique?,return first max
{HBTInt tmp,i;
	tmp=0;
	for(i=0;i<Len;i++)
	{
		if(vec[tmp]<vec[i])tmp=i;
	}
	return tmp;//max ID
}
HBTInt max_of_vec(HBTReal *vec,HBTInt Len)//what if the max is not unique?,return first max
{HBTInt tmp,i;
	tmp=0;
	for(i=0;i<Len;i++)
	{
		if(vec[tmp]<vec[i])tmp=i;
	}
	return tmp;//max ID
}

HBTInt Dmin_of_vec(double *vec,HBTInt Len)//what if the max is not unique???????
{HBTInt tmp,i;
	tmp=0;
	for(i=0;i<Len;i++)
	{
		if(vec[tmp]>vec[i])tmp=i;
	}
	return tmp;//min ID
}
HBTInt Fmin_of_vec(float *vec,HBTInt Len)//what if the max is not unique?,return first min
{HBTInt tmp,i;
	tmp=0;
	for(i=0;i<Len;i++)
	{
		if(vec[tmp]>vec[i])tmp=i;
	}
	return tmp;//min ID
}
HBTInt min_of_vec(HBTReal *vec,HBTInt Len)//what if the max is not unique?,return first min
{HBTInt tmp,i;
	tmp=0;
	for(i=0;i<Len;i++)
	{
		if(vec[tmp]>vec[i])tmp=i;
	}
	return tmp;//min ID
}

HBTInt insert_to_array(HBTInt a, HBTInt len, HBTInt *ascd_arr)//insert into an ascending arr
{
	HBTInt i;
	//~ if (a<ascd_arr[0])
	//~ {
		//~ fprintf(logfile,"error! too small to insert %d\n",a);
		//~ exit(1);
	//~ }
	for(i=0;i<len;i++)
	{ 
		if(ascd_arr[i]>a)
		return i-1;
	}
	return len-1;
}	

HBTInt prepare_ind2halo(CATALOGUE *A)
{HBTInt i,haloid,pid;
//	A->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);

	#pragma omp parallel 
		{
	#pragma omp for
	for(i=0;i<NP_DM;i++)//initialization
	{
		A->ID2Halo[i]=-1;/*return -1 if the PID does not belong to a Halo,
								i.e,we consider the backgroud as a halo with haloid=-1; 
								note that this only make sense when we try to find host for bound structures */
	}
	#pragma omp for private(i,pid,haloid)
	for(haloid=0;haloid<A->Ngroups;haloid++)
	{
		for(i=0;i<A->Len[haloid];i++)
		{
			pid=A->PIDorIndex[A->Offset[haloid]+i];//Pindex ranges [0,NP_DM);
			A->ID2Halo[pid]=haloid;//haloIDs begins from id=0
		}
	}
		}
	return 1; // means id2halo ready;
}

HBTInt * prepare_sub2halo(SUBCATALOGUE *SubCat)
{
	HBTInt i,haloid,subid;
	HBTInt *table_s2h;
	table_s2h=mymalloc(sizeof(HBTInt)*SubCat->Nsubs);	
	//prepare sub2halo table
	#pragma omp parallel 
	{
	#pragma omp for nowait private(haloid,i) //schedule(dynamic)
	for(haloid=0;haloid<SubCat->Ngroups;haloid++)
	{
		for(i=0;i<SubCat->GrpLen_Sub[haloid];i++)
		{
			table_s2h[SubCat->GrpOffset_Sub[haloid]+i]=haloid;
		}
	}
	//don't forget the quasihalos
	#pragma omp for
	for(subid=SubCat->Nsubs-SubCat->NQuasi;subid<SubCat->Nsubs;subid++)
		table_s2h[subid]=-1;//-1means background halo
	}
	return table_s2h;
}


HBTInt check_dup(HBTInt l1, HBTInt *p1,HBTInt l2, HBTInt *p2)
{
	HBTInt i,j,n;
	n=0;
	for(i=0;i<l1;i++)
		for(j=0;j<l2;j++)
		{
			if(p1[i]==p2[j])
				n++;
		}
	return n;
}

HBTReal position_modulus(HBTReal x)
{//shift the positions to within [0,boxsize)
	HBTReal y;
	if(x>=0&&x<BOXSIZE) return x;
	y=x/BOXSIZE;
	return (y-floor(y))*BOXSIZE;
}
	
HBTReal distance(HBTReal x[3],HBTReal y[3])
{
	HBTReal dx[3];
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	#ifdef PERIODIC_BDR
	dx[0]=NEAREST(dx[0]);
	dx[1]=NEAREST(dx[1]);
	dx[2]=NEAREST(dx[2]);
	#endif
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}

void center_of_mass(HBTReal CoM[3], HBTInt *PInd, HBTInt np, HBTReal PPos[][3])
{//center of mass function which properly handles periodic boundary condition
	HBTInt i,j;
	double sx[3],origin[3];
	
	if(0==np) return;
	
	sx[0]=sx[1]=sx[2]=0.;
	#ifdef PERIODIC_BDR
	for(j=0;j<3;j++)
		origin[j]=PPos[PInd[0]][j];
	#endif		
	for(i=0;i<np;i++)
	for(j=0;j<3;j++)
	#ifdef PERIODIC_BDR
		sx[j]+=NEAREST(PPos[PInd[i]][j]-origin[j]);
	#else
		sx[j]+=Pdat.Pos[PInd[i]][j];
	#endif
	for(j=0;j<3;j++)
	{
		sx[j]/=np;
		#ifdef PERIODIC_BDR
		sx[j]+=origin[j];
		#endif
		CoM[j]=sx[j];
	}
}
 
HBTReal comoving_virial_radius(HBTInt mass)
{	HBTReal Hratio,scaleF,virialF,x,OmegaZ;
	scaleF=header.time;
	Hratio=header.Hz/HUBBLE0;
	#ifdef OMEGA0
	OmegaZ=OMEGA0/(scaleF*scaleF*scaleF)/Hratio/Hratio;
	#else
	OmegaZ=header.Omega0/(scaleF*scaleF*scaleF)/Hratio/Hratio;
	#endif
	x=OmegaZ-1;
	virialF=18.0*3.1416*3.1416+82.0*x-39.0*x*x;//<Rho_vir>/Rho_cri
	return pow(2.0*G*mass*header.mass[1]/virialF/header.Hz/header.Hz,1.0/3)/scaleF;
}
void *mymalloc(size_t n)
{void * mem;
	if(n)
	{
		if(!(mem=malloc(n)))
		{fprintf(logfile,"failed to allocate memory for %zd bytes.\n", n);fflush(logfile);
		exit(1);
		}
	}
	else
	{
		mem=NULL;
	}
	return mem;
}
void myfree(void *mem)
{
	if(mem!=NULL)
	free(mem);
}


#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((unsigned) (x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))
//bit shift operation is invalid for 8byte+ data

void swap_Nbyte(void *data2swap,size_t nel,size_t mbyte)
/*This function is used to switch endian, for data2swap[nel] with element size mbyte*/
{
  size_t i,j;
  char *data, *old_data;//by definition, sizeof(char)=1, one byte
  
  data=(char *)data2swap;
  
  switch(mbyte)
  {
	  case 1 :break;
	  case 2 :
	  		for(j=0;j<nel;j++)
			FIX_SHORT(data[j*2]);
			break;
	  case 4 :
	  		for(j=0;j<nel;j++)
			FIX_LONG(data[j*4]);
			break;
	  default :
			old_data=mymalloc(mbyte);
			for(j=0;j<nel;j++)
			{
			  memcpy(&old_data[0],&data[j*mbyte],mbyte);
			  for(i=0;i<mbyte;i++)
				{
				  data[j*mbyte+i]=old_data[mbyte-i-1];
				}
			}
			myfree(old_data);
	}
}
size_t fread_swap(void *buf,size_t Nsize,size_t Nbuf,FILE *fp, int FlagByteSwap)
{
	size_t Nread;
	Nread=fread(buf,Nsize,Nbuf,fp);
	if(FlagByteSwap)
	swap_Nbyte(buf,Nbuf,Nsize);
	return Nread;
}

void spherical_basisD(double dx[3],double er[3],double et[3],double ef[3])
/*er: radial;
 *et: azimuthal,(0~2*pi)
 *ef: elevation (0~pi)
 * */ 
{
	int i;
	double dr,dxy2,modet,modef;//unit vectors of spherical coord.
	dxy2=dx[0]*dx[0]+dx[1]*dx[1];
	dr=sqrt(dxy2+dx[2]*dx[2]);
	for(i=0;i<3;i++) er[i]=dx[i]/dr;
	modet=sqrt(dxy2);
	if(0.==modet)
	{
		et[0]=1.;
		et[1]=et[2]=0.;
		ef[1]=1.;
		ef[0]=ef[2]=0.;
	}
	else
	{
		et[0]=-dx[1]/modet;
		et[1]=dx[0]/modet;
		et[2]=0;
		if(0.==dx[2])
		{
		ef[0]=ef[1]=0.;
		ef[2]=1.;	
		}
		else
		{
		modef=sqrt(dxy2+dxy2*dxy2/dx[2]/dx[2]);
		ef[0]=-dx[0]/modef;
		ef[1]=-dx[1]/modef;
		ef[2]=dxy2/dx[2]/modef;
		}
	}
}

HBTReal vec_prod(HBTReal *a,HBTReal *b,HBTInt dim)
{
	HBTInt i;
	HBTReal c;
	c=0.;
	for(i=0;i<dim;i++)
		c+=a[i]*b[i];
	return c;
}
void vec_cross(HBTReal a[3],HBTReal b[3],HBTReal c[3])
{
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[0]*b[2]-a[2]*b[0];
	c[2]=a[0]*b[1]-a[1]*b[0];
}
HBTInt try_readfile(char * filename)
{// return 1 when file can be read; 0 otherwise.
	FILE *fp;
	if(fp=fopen(filename,"r"))
	{
		 fclose(fp);
		 return 1;
	 }
	return 0;
}
HBTInt count_lines(char * filename)
{// return number of lines
	FILE *fp;
	if(!(fp=fopen(filename,"r")))
	{
		 return -1;
	 }
	 HBTInt i=0;
	 char ch;
	 while(1)
	 {
		 ch=fgetc(fp);
		 if(feof(fp))
		 {
			fseek(fp,-sizeof(char),SEEK_END);
			ch=fgetc(fp);
			if(ch!='\n') i++; //in case the last line is not terminated with '\n'
			break;
		 }
		 if(ch=='\n')
			i++;
	 }
	 fclose(fp);
	return i;
}
/* From Numerical Recipes*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
HBTReal psort(HBTInt k, HBTInt numel, HBTReal arr[])
{/* partition sort: sort arr into two parts: k small elements and n-k large elements,
  * with arr[k-1] as the k-th smallest value , also the return value.  
  * slightly modified from the original version to have arr[0~n-1] rather than arr[1~n]*/
	HBTInt i,ir,j,l,mid;
	HBTReal a,temp;

	arr--;//so that it's accessed with arr[1~n]
	l=1;
	ir=numel;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef SWAP

/*-------linkedlists-------------*/
HBTReal * GetArrPos(HBTInt pid,void *data)
{//pass Pdat.Pos to data
	return ((HBTxyz *) data)[pid];
}
HBTReal * GetSubCoM(HBTInt subid,void *data)
{//pass SubCat.Property to data
	return ((struct SubProperty *) data)[subid].CoM;
}
HBTInt linklist_round_gridid(HBTInt i,HBTInt ndiv)
{//to correct for rounding error near boundary 
	if(i<0) return 0;
	if(i>=ndiv) return ndiv-1;
	return i;
}
HBTInt linklist_shift_gridid(HBTInt i,HBTInt ndiv)
{//to correct for periodic conditions; 
//only applicable when def PERIODIC_BDR and ll.UseFullBox=1
	i=i%ndiv;
	if(i<0) i+=ndiv;
	return i;
}
HBTInt linklist_fix_gridid(HBTInt i, LINKLIST *ll)
{
	#ifdef PERIODIC_BDR
	if(ll->UseFullBox)
	return linklist_shift_gridid(i,ll->ndiv);	
	#endif
	return linklist_round_gridid(i,ll->ndiv);
}
HBTInt linklist_get_hoc(LINKLIST *ll, HBTInt i,HBTInt j,HBTInt k)
{
	return ll->hoc[i+j*ll->ndiv+k*ll->ndiv*ll->ndiv];
}
HBTInt linklist_get_hoc_safe(LINKLIST *ll, HBTInt i,HBTInt j,HBTInt k)
{//force fixing of gridids
	#define FIXGRID(i) linklist_fix_gridid(i,ll)
	return ll->hoc[FIXGRID(i)+FIXGRID(j)*ll->ndiv+FIXGRID(k)*ll->ndiv*ll->ndiv];
}

void make_linklist(LINKLIST *ll, HBTInt np,HBTInt ndiv, void *PosData, 
								AccessPosFunc *GetPos, HBTInt UseFullBox)
{
	HBTInt i,j,grid[3];
	HBTInt ind,ndiv2;
	HBTReal x;
	//~ float range[3][2],step[3];
	printf("creating linked list..\n");
	
	ndiv2=ndiv*ndiv;
	ll->ndiv=ndiv;
	ll->np=np;
	ll->UseFullBox=UseFullBox;
	ll->hoc=mymalloc(sizeof(HBTInt)*ndiv*ndiv*ndiv);
	ll->list=mymalloc(sizeof(HBTInt)*np);
	ll->PosData=PosData;
	ll->GetPos=GetPos;
	/*determining enclosing cube*/
	if(UseFullBox)
	{
		for(i=0;i<3;i++)
		{
			ll->range[i][0]=0.;
			ll->range[i][1]=BOXSIZE;
		}
		for(j=0;j<3;j++)
			ll->step[j]=BOXSIZE/ndiv;	
	}
	else
	{
		for(i=0;i<3;i++)
			for(j=0;j<2;j++)
				ll->range[i][j]=GetPos(0,PosData)[i];
		for(i=1;i<np;i++)
			for(j=0;j<3;j++)
			{
				x=GetPos(i,PosData)[j];
				if(x<ll->range[j][0])
					ll->range[j][0]=x;
				else if(x>ll->range[j][1])
					ll->range[j][1]=x;
			}
		for(j=0;j<3;j++)
			ll->step[j]=(ll->range[j][1]-ll->range[j][0])/ll->ndiv;
	}
	/*initialize hoc*/
	HBTInt *phoc=ll->hoc;
	for(i=0;i<ndiv*ndiv*ndiv;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((GetPos(i,PosData)[j]-ll->range[j][0])/ll->step[j]);
			grid[j]=linklist_fix_gridid(grid[j],ll);
		}
		ind=grid[0]+grid[1]*ndiv+grid[2]*ndiv2;
		ll->list[i]=ll->hoc[ind];
		ll->hoc[ind]=i; /*use hoc[ind] as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
}

void free_linklist(LINKLIST *ll)
{
	myfree(ll->hoc);
	myfree(ll->list);
	ll->ndiv=0;
	ll->np=0;
	ll->PosData=NULL;
	ll->GetPos=NULL;
}

HBTInt * linklist_search_sphere(LINKLIST *ll, HBTReal radius, HBTReal searchcenter[3], HBTInt *N_max_and_found)
{
	HBTReal dr;
	HBTInt i,j,k,pid,*PIDfound,nfound,nmax;
	int subbox_grid[3][2];
	
	nmax=*N_max_and_found;
	PIDfound=mymalloc(sizeof(HBTInt)*nmax);
	nfound=0;
		
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((searchcenter[i]-radius-ll->range[i][0])/ll->step[i]);
	subbox_grid[i][1]=floor((searchcenter[i]+radius-ll->range[i][0])/ll->step[i]);
#ifndef PERIODIC_BDR //do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
	subbox_grid[i][0]=linklist_fix_gridid(subbox_grid[i][0],ll);
	subbox_grid[i][1]=linklist_fix_gridid(subbox_grid[i][1],ll);
#endif	
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc_safe(ll,i,j,k); //in case the grid-id is out of box, in the periodic case
// 				pid=linklist_get_hoc(ll,i,j,k);
				while(pid>=0)
				{
					dr=distance(ll->GetPos(pid,ll->PosData),searchcenter);
					if(dr<radius)
					{
						if(nfound==nmax)
						{
							nmax*=2;
							PIDfound=realloc(PIDfound,sizeof(HBTInt)*nmax);
						}
						PIDfound[nfound]=pid;					
						nfound++;
					}
					pid=ll->list[pid];
				}
			}
	*N_max_and_found=nfound;		
	return PIDfound;		
}



void av_center(int *Arr, int Np, HBTReal CoM[3])
{
	int i,j;
	double com[3]={0.};
	for(i=0;i<Np;i++)
	{
		for(j=0;j<3;j++)
			com[j]+=Pdat.Pos[Arr[i]][j];
	}
	for(j=0;j<3;j++)
			CoM[j]=com[j]/(double)Np;
}

#define RCONTRACT 0.8
#define RTOLERATE 4
#define NCoMMin 10

void moving_center(int *Arr, int Np, HBTReal CoM[3])
{
	int i,j,Ncore,*ArrNew,*ArrOld;
	HBTReal *r, rmax, Cen[3];
	
	ArrOld=mymalloc(sizeof(int)*Np);
	ArrNew=mymalloc(sizeof(int)*Np);
	r=mymalloc(sizeof(HBTReal)*Np);
	rmax=0;
	
	/*clean the data for low-resolution particles, meant for AHF*/
	Ncore=0;
	for(i=0;i<Np;i++)
	{
		if(Arr[i]>=0&&Arr[i]<NP_DM)
		{
			ArrOld[Ncore]=Arr[i];
			Ncore++;
		}
	}
	
	av_center(ArrOld,Ncore,Cen); //init Cen

	for(i=0;i<Ncore;i++)//init rmax
	{
		r[i]=distance(Pdat.Pos[ArrOld[i]], Cen);
		if(r[i]>rmax) rmax=r[i];
	}
	rmax*=RCONTRACT; //contract
	for(i=0,j=0;i<Ncore;i++) 
	{
		if(r[i]<rmax)
		{
			ArrNew[j]=ArrOld[i];
			j++;
		}
	}
	free(r);
	Ncore=j; //contracted
	av_center(ArrNew,Ncore,CoM); //new CoM
	while(distance(CoM,Cen)>RTOLERATE*SofteningHalo)
	{
		for(i=0;i<Ncore;i++)
			ArrOld[i]=ArrNew[i];
		for(i=0;i<3;i++)
			Cen[i]=CoM[i];
		rmax*=RCONTRACT;
		for(i=0,j=0;i<Ncore;i++)
		{
			if(distance(Pdat.Pos[ArrOld[i]], Cen)<rmax)
			{
				ArrNew[j]=ArrOld[i];
				j++;
			}
		}
		Ncore=j;
		if(Ncore<NCoMMin)
		{
			printf("warning: CoM has not converged for mass=%d\n", Np);
			break;
		}
		av_center(ArrNew,Ncore,CoM);
	}	
	free(ArrNew);
	free(ArrOld);
}
