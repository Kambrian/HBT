#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void multiprint()
{
	int i;
	for(i=0;i<5;i++)
	printf("from thread %d, i=%d\n",omp_get_thread_num(),i);
}
void multiprint2()
{
	int i;
	for(i=0;i>-5;i--)
	printf("from thread %d, i=%d\n",omp_get_thread_num(),i);
}	
int main()
{
	#pragma omp parallel sections num_threads(5)
	{
		#pragma omp section
		multiprint();
		#pragma omp section
		multiprint2();
		#pragma omp section
		printf("from thread %d, single\n",omp_get_thread_num());
	}
	return 0;
}
