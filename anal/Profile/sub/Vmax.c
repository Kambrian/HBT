/* to extract maximum circular velocity Vmax and comoving Rmax for subhalos *
 * also gives comoving half-mass radius rhalf *
 * 1,2,3*sigma percentile radius rsig,r2sig,r3sig
 * poisson percentile radius rpoisson ((N-sqrt(N))/N percentile)
 * output: RmaxVmax_%d.CoM or RmaxVmax_%d.MBD, using CoM or MostBound as center *
 *         CenVCen_%d.***, Center (comoving) and VCenters (physical), 
 * 	   can also output Mvir and Rvir using only particles within the subhalo
 * units: standard units */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define CEN_COM 0
#define CEN_MBD 1
#define CEN_MPT 2  //when using this, need to define HALO_PARA in the parameter file to make the tree-code thread-safe

#define CEN_TYPE CEN_COM
#define OUTPUT_SUBVIR

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

static int compare_real(const void *a, const void *b)
{
  if((*(HBTReal *) a) > (*(HBTReal *) b))
    return +1;

  if((*(HBTReal *) a) < (*(HBTReal *) b))
    return -1;

  return 0;
}
HBTInt min_pot_particle(HBTInt subid);
void halo_virial_factor(HBTReal virialF[3]);
void halo_virial_radius(HBTReal Mvir[3], HBTReal Rvir[3], HBTReal r_sorted[], HBTInt np, HBTReal virialF[3]);

SUBCATALOGUE SubCat;

int main()
{
char outputdir[1024];
	
CATALOGUE Cat;
HBTInt Nsnap,i,subid,sublen,imax,iminpot;
HBTReal *cen;
HBTReal *r,*rmax,*vmax,*rhalf,*rsig,*r2sig,*r3sig,*rpoisson, (*Mvir)[3], (*Rvir)[3],virialF[3];  
HBTReal *v,unit,sqa;
HBTxyz *centers, *vcenters;

char buf[1024];
FILE *fp;
logfile=stdout;
sprintf(outputdir,"%s/profile/",SUBCAT_DIR);	
mkdir(outputdir,0755);

load_particle_header(0,SNAPSHOT_DIR);
unit=sqrt(G*header.mass[1]);

for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&SubCat,-2);
	free_PIDHash();
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(header.time);
	#endif
		
	rmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	vmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	rhalf=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	rsig=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	r2sig=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	r3sig=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	rpoisson=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	centers=mymalloc(sizeof(HBTxyz)*SubCat.Nsubs);
	vcenters=mymalloc(sizeof(HBTxyz)*SubCat.Nsubs);
	Mvir=mymalloc(sizeof(HBTReal)*SubCat.Nsubs*3);
	Rvir=mymalloc(sizeof(HBTReal)*SubCat.Nsubs*3);
#ifdef OUTPUT_SUBVIR
	halo_virial_factor(virialF);
#endif
	#pragma omp parallel for private(subid,cen,sublen,r,v,i,imax,iminpot)	
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		if(SubCat.SubLen[subid])
		{
			#if CEN_TYPE==CEN_COM
			cen=SubCat.Property[subid].CoM;
			for(i=0;i<3;i++) vcenters[subid][i]=SubCat.Property[subid].VCoM[i];
			#elif CEN_TYPE==CEN_MBD
			cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
			for(i=0;i<3;i++) vcenters[subid][i]=Pdat.Vel[SubCat.PSubArr[subid][0]][i]*sqa;
			#elif CEN_TYPE==CEN_MPT
			iminpot=min_pot_particle(subid);
			cen=Pdat.Pos[iminpot];
			for(i=0;i<3;i++) vcenters[subid][i]=Pdat.Vel[iminpot][i]*sqa;
			#endif
			for(i=0;i<3;i++) centers[subid][i]=cen[i];
			sublen=SubCat.SubLen[subid];
			r=mymalloc(sizeof(HBTReal)*sublen);	
			v=mymalloc(sizeof(HBTReal)*sublen);
			for(i=0;i<sublen;i++)
			   r[i]=distance(cen,Pdat.Pos[SubCat.PSubArr[subid][i]]);	
			qsort(r,sublen,sizeof(HBTReal),compare_real);
			#ifdef OUTPUT_SUBVIR
			halo_virial_radius(Mvir[subid],Rvir[subid],r,sublen,virialF);
			#endif
			for(i=0;i<sublen;i++)
			{
				if(r[i]<SofteningHalo) r[i]=SofteningHalo; //resolution
				v[i]=sqrt((HBTReal)(i+1)/r[i]);
			}
			imax=max_of_vec(v,sublen);
			rmax[subid]=r[imax];
			vmax[subid]=v[imax]*unit;
			rhalf[subid]=r[sublen/2];
			rsig[subid]=r[(HBTInt)(sublen*0.683)];
			r2sig[subid]=r[(HBTInt)(sublen*0.955)];
			r3sig[subid]=r[(HBTInt)(sublen*0.997)];
			rpoisson[subid]=r[sublen-(HBTInt)sqrtf(sublen)];
			myfree(r);
			myfree(v);
		}
		else
		{
			rmax[subid]=0.;
			vmax[subid]=0.;
			rhalf[subid]=0.;
			rsig[subid]=0.;
			r2sig[subid]=0.;
			r3sig[subid]=0.;
			rpoisson[subid]=0.;
			for(i=0;i<3;i++)
			{
			centers[subid][i]=0.;
			vcenters[subid][i]=0.;
			Mvir[subid][i]=0.;
			Rvir[subid][i]=0.;
			}
		}
	}

#if CEN_TYPE==CEN_COM
	sprintf(buf,"%s/RmaxVmax_"HBTIFMT".COM",outputdir,Nsnap);
#elif CEN_TYPE==CEN_MBD
	sprintf(buf,"%s/RmaxVmax_"HBTIFMT".MBD",outputdir,Nsnap);
#elif CEN_TYPE==CEN_MPT
	sprintf(buf,"%s/RmaxVmax_"HBTIFMT".MPT",outputdir,Nsnap);
#endif	
	myfopen(fp,buf,"w");	
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);
	fwrite(rmax,sizeof(HBTReal),SubCat.Nsubs,fp);
	fwrite(vmax,sizeof(HBTReal),SubCat.Nsubs,fp);
	fwrite(rhalf,sizeof(HBTReal),SubCat.Nsubs,fp);
	fwrite(rsig,sizeof(HBTReal),SubCat.Nsubs,fp);
	fwrite(r2sig,sizeof(HBTReal),SubCat.Nsubs,fp);
	fwrite(r3sig,sizeof(HBTReal),SubCat.Nsubs,fp);
	fwrite(rpoisson,sizeof(HBTReal),SubCat.Nsubs,fp);
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);
	fclose(fp);
	
#if CEN_TYPE==CEN_COM
	sprintf(buf,"%s/CenVCen_"HBTIFMT".COM",outputdir,Nsnap);
#elif CEN_TYPE==CEN_MBD
	sprintf(buf,"%s/CenVCen_"HBTIFMT".MBD",outputdir,Nsnap);
#elif CEN_TYPE==CEN_MPT
	sprintf(buf,"%s/CenVCen_"HBTIFMT".MPT",outputdir,Nsnap);
#endif
	myfopen(fp,buf,"w");	
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);
	for(i=0;i<SubCat.Nsubs;i++)
	fwrite(centers[i],sizeof(HBTReal),3,fp);
	for(i=0;i<SubCat.Nsubs;i++)
	fwrite(vcenters[i],sizeof(HBTReal),3,fp);
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);
	fclose(fp);

#ifdef OUTPUT_SUBVIR	
#if CEN_TYPE==CEN_COM
	sprintf(buf,"%s/SubMvirRvir_"HBTIFMT".COM",outputdir,Nsnap);
#elif CEN_TYPE==CEN_MBD
	sprintf(buf,"%s/SubMvirRvir_"HBTIFMT".MBD",outputdir,Nsnap);
#elif CEN_TYPE==CEN_MPT
	sprintf(buf,"%s/SubMvirRvir_"HBTIFMT".MPT",outputdir,Nsnap);
#endif
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);
	for(i=0;i<SubCat.Nsubs;i++)
	fwrite(Mvir[i],sizeof(HBTReal),3,fp);
	for(i=0;i<SubCat.Nsubs;i++)
	fwrite(Rvir[i],sizeof(HBTReal),3,fp);
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);
	fclose(fp);
#endif
	
	myfree(rmax);
	myfree(vmax);
	myfree(rhalf);
	myfree(rsig);
	myfree(r2sig);
	myfree(r3sig);
	myfree(rpoisson);
	myfree(centers);
	myfree(vcenters);
	myfree(Mvir);
	myfree(Rvir);
	erase_sub_catalogue(&SubCat);
	free_particle_data();
}

return 0;
}

HBTInt load_CenVCen(HBTxyz *Cen,HBTxyz *VCen, HBTInt Nsnap)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy,i;
	
	if(!Cen||!VCen)
	{
		printf("error: allocate cen, vcen first \n");
		exit(1);
	}
	#if CEN_TYPE==CEN_COM
	sprintf(buf,"%s/profile/CenVCen_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MBD
	sprintf(buf,"%s/profile/CenVCen_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MPT
	sprintf(buf,"%s/profile/CenVCen_"HBTIFMT".MPT",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	for(i=0;i<Nsubs;i++)
	fread(Cen[i],sizeof(HBTReal),3,fp);
	for(i=0;i<Nsubs;i++)
	fread(VCen[i],sizeof(HBTReal),3,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
	return Nsubs;
}

HBTInt load_RmaxVmax(HBTReal *rmax,HBTReal *vmax, HBTReal *rhalf, HBTInt Nsnap)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	//~ HBTReal *rsig, *r2sig, *r3sig, *rpoisson; /* declare this as input var if you want them!!! */
	
	if(!rmax||!vmax||!rhalf)
	{
		printf("error: allocate rmax , vmax and rhalf first \n");
		exit(1);
	}
	#if CEN_TYPE==CEN_COM
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MBD
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MPT
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MPT",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	fread(rmax,sizeof(HBTReal),Nsubs,fp);
	fread(vmax,sizeof(HBTReal),Nsubs,fp);
	fread(rhalf,sizeof(HBTReal),Nsubs,fp);
	fseek(fp,sizeof(HBTReal)*Nsubs*4,SEEK_CUR);
	//~ fread(rsig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r2sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r3sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(rpoisson,sizeof(HBTReal),Nsubs,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
	return Nsubs;
}

HBTInt min_pot_particle(HBTInt subid)
{
	double potmin,pot;
	HBTInt ipotmin,i;
	
		if(0==SubCat.SubLen[subid]) 
				return -1;
		
		potmin=0.;ipotmin=0;
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubCat.SubLen[subid],SubCat.SubLen[subid]);
		maketree(SubCat.SubLen[subid],SubCat.PSubArr[subid],Pdat.Pos);
		for(i=0;i<SubCat.SubLen[subid];i++)
		{
			pot=tree_treeevaluate_potential(Pdat.Pos[SubCat.PSubArr[subid][i]],SubCat.PSubArr[subid],Pdat.Pos);
			if(potmin>pot)
			{
			ipotmin=i;
			potmin=pot;
			}
		}
		tree_tree_free();
		return SubCat.PSubArr[subid][ipotmin];
}

void halo_virial_factor(HBTReal virialF[3])
{
	HBTReal Hratio,scaleF,virialF_tophat,virialF_c200,virialF_b200,x,OmegaZ;
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
void halo_virial_radius(HBTReal Mvir[3], HBTReal Rvir[3], HBTReal r_sorted[], HBTInt np, HBTReal virialF[3])
{//works on sorted radius list r_sorted[]
  HBTReal tol=1e-5;
  HBTInt i,ndiv;
  HBTReal *r,rvir,rdiv;
  HBTInt virtype;
  
  for(virtype=0;virtype<3;virtype++)
  {
    ndiv=np;//guess mass
    r=r_sorted;
    rdiv=r[ndiv];
    rvir=pow(2.0*G*ndiv*header.mass[1]/virialF[virtype]/header.Hz/header.Hz,1.0/3)/header.time;//guess radius
    while(1)
    {
      if(rdiv>rvir)//reduce mass guess
      {
	for(i=ndiv-1;i>=0;i--)
	{
	  if(r[i]<rvir) break;
	}
	ndiv=i+1;
      }
      else if(rdiv<rvir) //increase mass guess
      {
	for(i=ndiv;i<np;i++)
	{
	  if(r[i]>=rvir) break;
	}
	ndiv=i;
      }
      
      rdiv=rvir;
      rvir=pow(2.0*G*ndiv*header.mass[1]/virialF[virtype]/header.Hz/header.Hz,1.0/3)/header.time;//recalibrate radius
      
      if(0==ndiv||ndiv==np||fabs((rvir-rdiv)/rvir)<tol) break; //converged
    }

    Rvir[virtype]=rvir;
    Mvir[virtype]=ndiv*header.mass[1];
  }
}