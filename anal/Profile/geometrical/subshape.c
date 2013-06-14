#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define GRPID_MAX 100
#define SUBLEN_MIN 100

//~ #define SUBFIND_DIR "/home/kambrain/data/8213/subcatS"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat;
	int i,subid,pid;
	float *cen,dx,dy,dz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,(*Q)[3];
	double *data;
	
	
	gsl_matrix_view m ;
	gsl_eigen_symm_workspace * w;
	gsl_vector *eval;


	FILE *fp;
	char buf[1024];
	int Nsnap;

	logfile=stdout;
	
	Nsnap=atoi(argv[1]);
	
	sprintf(buf,"%s/anal/sub_shape_%d",SUBCAT_DIR,Nsnap);
	if(!(fp=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fresh_ID2Index(&SubCat,-2);
	#ifdef GRPID_MAX
	SubCat.Nsubs=SubCat.GrpOffset_Sub[GRPID_MAX];
	#endif
	Q=mymalloc(sizeof(float)*3*SubCat.Nsubs);
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
		//~ cen=SubCat.Property[subid].CoM;
		cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
		for(i=0;i<SubCat.SubLen[subid];i++)
		{
			pid=SubCat.PSubArr[subid][i];
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
		free(data);
		gsl_vector_free (eval);
		gsl_eigen_symm_free (w);
	}
	for(subid=0;subid<SubCat.Nsubs;subid++)
	fprintf(fp,"%g,%g,%g\n",Q[subid][0],Q[subid][1],Q[subid][2]);
	fclose(fp);
	return 0;
}
