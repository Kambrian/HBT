#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int main()
{
	float Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Q[3];
	double *data;
	int i;
	
	gsl_matrix_view m ;
	gsl_eigen_symm_workspace * w;
	gsl_vector *eval;


		Ixx=Iyy=Izz=Ixy=Ixz=Iyz=100.;
		data=malloc(sizeof(double)*9);
		data[0]=Ixx;data[1]=Ixy;data[2]=Ixz;data[3]=Ixy;data[4]=Iyy;data[5]=Iyz;data[6]=Ixz;data[7]=Iyz;data[8]=Izz;
		m= gsl_matrix_view_array (data, 3, 3);
		w = gsl_eigen_symm_alloc (3);
		eval = gsl_vector_alloc (3);
		gsl_eigen_symm(&m.matrix, eval, w);
		for(i=0;i<3;i++)
			Q[i]=gsl_vector_get(eval,i);
		free(data);
		gsl_vector_free (eval);
		gsl_eigen_symm_free (w);
	printf("%g,%g,%g\n",Q[0],Q[1],Q[2]);
	return 0;
}
