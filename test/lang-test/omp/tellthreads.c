//to tell how many threads are in use
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char** argv)
{
	//~ omp_set_dynamic(1);
	printf("total number of threads = %d\n", omp_get_num_threads());
	printf("omp dynamic thread adjustion is %d\n",omp_get_dynamic());
	#pragma omp parallel
	{
		#pragma omp single
		{
		printf("total number of threads = %d\n", omp_get_num_threads());
		printf("#####################\n");
		}
		printf("thread %d out of %d\n",omp_get_thread_num(),omp_get_num_threads());
	}
	printf("#####################\n");
	return 0;
}
	
