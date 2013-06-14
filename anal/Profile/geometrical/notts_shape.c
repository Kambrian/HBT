#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define SUBLEN_MIN 10

typedef struct
{
	int Nsubs;
	int Nids;
	int *SubLen;
	int *SubOffset;
	int *SubName;
	int *SubArr;
	HBTxyz * CoM;
} NOTTSCATALOGUE;

SUBCATALOGUE SCat;
NOTTSCATALOGUE SubCat;
int load_notts_catalog(char * code);
void av_center(int *Arr, int Np, HBTReal CoM[3]);
void moving_center(int *Arr, int Np, HBTReal CoM[3]);
static int compare_dbl(const void *a, const void *b)//used to sort desendent id in descending order
{
  if((* (double *) a) > (*(double *) b))
    return -1;

  if((*(double *) a) < (*(double *) b))
    return +1;

  return 0;
}
int main(int argc,char **argv)
{
	int i,subid,pid;
	double dx,dy,dz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,(*Q)[3];
	double *data;
	HBTReal CenHalo[3]={57060.4, 52618.6, 48704.8},*d,*cen;
	
	gsl_matrix_view m ;
	gsl_eigen_symm_workspace * w;
	gsl_vector *eval;


	FILE *fp;
	char buf[1024],*code;
	int Nsnap;

	code=argv[1];
	logfile=stdout;
	Nsnap=MaxSnap-1;
	#ifdef PERIODIC_BDR
	printf("PERIODIC!\n");
	#endif
	
	sprintf(buf,"/home/sgnworkshop/DATA/results/shape/shape_raw_AqA4_%s",code);
	
	myfopen(fp,buf,"w");
	fprintf(fp,"I1,I2,I3,r,m,id\n");
	load_notts_catalog(code);
//    load_sub_catalogue(Nsnap,&SCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(SubCat.SubArr,SubCat.Nids);
//	fresh_ID2Index(&SCat,FRSH_SUBCAT);
	free_PIDHash();
	
/*	HBTReal CoM0[3], CoM1[3];
	subid=0;
	av_center(SCat.PSubArr[subid],SCat.SubLen[subid],CoM0);
	moving_center(SCat.PSubArr[subid],SCat.SubLen[subid],CoM1);
	printf("%g,%g,%g: %g\n", SCat.Property[subid].CoM[0],SCat.Property[subid].CoM[1],SCat.Property[subid].CoM[2],distance(SCat.Property[subid].CoM,SCat.Property[subid].CoM));
	printf("%g,%g,%g: %g\n", CoM0[0],CoM0[1],CoM0[2],distance(CoM0,SCat.Property[subid].CoM));
	printf("%g,%g,%g: %g\n", CoM1[0],CoM1[1],CoM1[2],distance(CoM1,SCat.Property[subid].CoM));
*/
		
	SubCat.CoM=mymalloc(sizeof(HBTReal)*3*SubCat.Nsubs);
	d=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		moving_center(SubCat.SubArr+SubCat.SubOffset[subid],SubCat.SubLen[subid],SubCat.CoM[subid]);
		d[subid]=distance(SubCat.CoM[subid],CenHalo);
	}
	
	Q=mymalloc(sizeof(double)*3*SubCat.Nsubs);
	#pragma omp parallel for private(subid,cen,i,pid,dx,dy,dz,data,m,w,eval,Ixx,Iyy,Izz,Ixy,Ixz,Iyz)
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		if(SubCat.SubLen[subid]<SUBLEN_MIN)
		{
			for(i=0;i<3;i++)
				Q[subid][i]=0;
			 continue;
		}
		Ixx=Iyy=Izz=Ixy=Ixz=Iyz=0.;
		cen=SubCat.CoM[subid];
		for(i=SubCat.SubOffset[subid];i<SubCat.SubOffset[subid]+SubCat.SubLen[subid];i++)
		{
			pid=SubCat.SubArr[i];
			if(pid<0||pid>=NP_DM) continue;
			dx=Pdat.Pos[pid][0]-cen[0];
			dy=Pdat.Pos[pid][1]-cen[1];
			dz=Pdat.Pos[pid][2]-cen[2];
			Ixx+=dx*dx;
			Iyy+=dy*dy;
			Izz+=dz*dz;
			Ixy+=dx*dy;
			Ixz+=dx*dz;
			Iyz+=dy*dz;
			//~ Ixx+=dy*dy+dz*dz;
			//~ Iyy+=dx*dx+dz*dz;
			//~ Izz+=dx*dx+dy*dy;
			//~ Ixy-=dx*dy;
			//~ Ixz-=dx*dz;
			//~ Iyz-=dy*dz;
		}
		data=mymalloc(sizeof(double)*9);
		data[0]=Ixx;data[1]=Ixy;data[2]=Ixz;data[3]=Ixy;data[4]=Iyy;data[5]=Iyz;data[6]=Ixz;data[7]=Iyz;data[8]=Izz;
		m= gsl_matrix_view_array (data, 3, 3);
		w = gsl_eigen_symm_alloc (3);
		eval = gsl_vector_alloc (3);
		gsl_eigen_symm(&m.matrix, eval, w);
		for(i=0;i<3;i++)
			Q[subid][i]=gsl_vector_get(eval,i);
		qsort(Q[subid],3,sizeof(double),compare_dbl);	
		free(data);
		gsl_vector_free (eval);
		gsl_eigen_symm_free (w);
	}
	for(subid=0;subid<SubCat.Nsubs;subid++)
	fprintf(fp,"%g,%g,%g,%g,%d,%d\n",Q[subid][0],Q[subid][1],Q[subid][2],d[subid],SubCat.SubLen[subid],SubCat.SubName[subid]);
	fclose(fp);
	return 0;
}

static int LINE_NUM=1;

int finish_read_line(FILE *fp)
{
	int status;
	status=1;
	fseek(fp,-1,SEEK_CUR); 
	if(fgetc(fp)!='\n') 
		status=seek_new_line(fp);
	else
		LINE_NUM++;
		
	return status;
}

int read_tag_line(int *subid, int *sublen, FILE *fp)
{
	int flag_end;
	flag_end=skip_comment(fp);
	if(flag_end<0)
	{
		printf("End of File.\n");
		return flag_end;
	}
	flag_end=get_new_int(subid,fp);
	if(flag_end<0) 
	{
		printf("error: no input on line %d\n",LINE_NUM);
		exit(1);	
	}
	flag_end=get_new_int(sublen,fp);
	if(flag_end<0) 
	{
		printf("error: SubLen expected on line %d\n",LINE_NUM);
		exit(1);	
	}
	flag_end=finish_read_line(fp);
	return flag_end;
}

int seek_new_line(FILE *fp)
{
	char tmp;
	tmp=fgetc(fp);
	while(tmp!='\n')
	{
		if(feof(fp))
			return -1;
		tmp=fgetc(fp);
	}
	LINE_NUM++;
	return 1;
}
int skip_comment(FILE *fp)
{
	//~ char tmp;
	while(!isdigit(fgetc(fp)))
	{
		if(feof(fp))
			return -1;
		seek_new_line(fp);
	}
	fseek(fp,-sizeof(char),SEEK_CUR);
	return 1;
}
int get_new_word(char *s,FILE *fp,int strlen_max)
{
	int n;
	n=0;
	while(!(isgraph(*s=fgetc(fp))))
	{
		if(s[0]=='\n')
		{
		LINE_NUM++;
		return -1;
		}
		if(feof(fp))
		return -2;
	}
	n++;
	s++;
	while(isgraph(*s=fgetc(fp)))
	{
		s++;
		n++;
		if(n==strlen_max-1)
		{
			*s='\0';
			return n;
		}
	}
	*s='\0';
	return n;
}

int get_new_int(int *i, FILE *fp)
{
	char s[10];
	int stat;
	if((stat=get_new_word(s,fp,10))>0)
	{
		(*i)=strtol(s,NULL,10);
		return 1;
	}
	
	return stat;	
}


int load_notts_catalog(char * code)
{
	char buf[1024];
	FILE *fp;
	int flag_end;
	int Nsub,Npart,NsubMax;
	
	sprintf(buf,"/home/sgnworkshop/jiaxin/data/%s-4",code);
	myfopen(fp,buf,"r");
	
	NsubMax=5000;
	SubCat.SubLen=mymalloc(sizeof(int)*NsubMax);
	SubCat.SubOffset=mymalloc(sizeof(int)*NsubMax);
	SubCat.SubName=mymalloc(sizeof(int)*NsubMax);
	SubCat.SubArr=mymalloc(sizeof(int)*NP_DM);
	SubCat.Nsubs=0;
	SubCat.Nids=0;
	Nsub=0;
	Npart=0;
	while(1)
	{
		if(Nsub==NsubMax)
		{
			NsubMax*=2;
			SubCat.SubLen=realloc(SubCat.SubLen,sizeof(int)*NsubMax);
			SubCat.SubOffset=realloc(SubCat.SubOffset,sizeof(int)*NsubMax);
			SubCat.SubName=realloc(SubCat.SubName,sizeof(int)*NsubMax);
		}	
		flag_end=read_tag_line(SubCat.SubName+Nsub,SubCat.SubLen+Nsub,fp);
		if(flag_end<0) break;
//		printf("%d,%d\n",SubCat.SubName[Nsub],SubCat.SubLen[Nsub]);
		SubCat.Nids+=SubCat.SubLen[Nsub];
		SubCat.SubOffset[Nsub]=Npart;
		while(Npart<SubCat.Nids)
		{
			skip_comment(fp);
			flag_end=get_new_int(SubCat.SubArr+Npart,fp);
//			printf("%d: %d\n",LINE_NUM, SubCat.SubArr[Npart]);
			if(flag_end<0)
			{
				printf("error: pid expected on line %d\n",LINE_NUM);
				exit(1);	
			}
			flag_end=finish_read_line(fp);
			Npart++;
		}
		Nsub++;
		if(flag_end<0) break;
	}
	SubCat.Nsubs=Nsub;
	SubCat.SubLen=realloc(SubCat.SubLen,sizeof(int)*Nsub);
	SubCat.SubOffset=realloc(SubCat.SubOffset,sizeof(int)*Nsub);
	SubCat.SubName=realloc(SubCat.SubName,sizeof(int)*Nsub);
	SubCat.SubArr=realloc(SubCat.SubArr,sizeof(int)*Npart);
	printf("finished reading %d lines\n", LINE_NUM);
	printf("%d subhalos with %d particles loaded\n",Nsub,Npart);
	return Npart;
}

void av_center(int *Arr, int Np, HBTReal CoM[3])
{
	int i,j;
	CoM[0]=CoM[1]=CoM[2]=0.;
	for(i=0;i<Np;i++)
	{
		for(j=0;j<3;j++)
			CoM[j]+=Pdat.Pos[Arr[i]][j];
	}
	for(j=0;j<3;j++)
			CoM[j]/=Np;
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

