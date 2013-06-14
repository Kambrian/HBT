//to create a long lasting process with give num of threads
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

int main(int argc, char** argv)
{
	int i,Nthreads;
	float j;
	if(argc!=2)
	{
		printf("usage: %s [Nthreads]\n", argv[0]);
		exit(1);
	}
	Nthreads=atoi(argv[1]);
	omp_set_num_threads(Nthreads);
	#pragma omp parallel private(i,j)
	while(1)
	{
		j=0;
		for(i=0;i<1024*1024;i++)
		j+=i*0.1;
	}
	return 0;
}
	
