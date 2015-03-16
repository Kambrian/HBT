#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main()
{
	int i,j;
	i=1;
	#pragma omp parallel
	{
	if(i==0)
	{
		#pragma omp for
		for(j=0;j<5;j++)
		printf("j=%d,hello from thread %d\n",j,omp_get_thread_num());
	}
	else
	#pragma omp single
	{
		printf("single hello\n");
	}
	}
return 0;
}
