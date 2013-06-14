//to test if dynamical adjustion have effect
//Result: dynamical adjustion is effective for icc v11.1.072
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

void task()
{
	int j;
	float k;
	k=0.;
	for(j=0;j<1024*1024;j++)
		k+=j*0.1+k*k-j*j*k;
}
int main(int argc, char** argv)
{
	int i,j,Nthreads;

	omp_set_dynamic(1);
	printf("omp dynamic thread adjustion is %d\n",omp_get_dynamic());
	for(i=0;i<100;i++)
	{
		#pragma omp parallel
		{
		#pragma omp single 
		printf("i=%d,Nthreads=%d,",i,omp_get_num_threads());	 
		task();
		#pragma omp single 
		printf(" %d\n",omp_get_num_threads());
		}
	}
	return 0;
}
	
