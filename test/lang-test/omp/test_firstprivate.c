#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
	
int main()
{
	int i, a;
	#pragma omp parallel num_threads(2)
	{
	#pragma omp single
	  for(i=0;i<10;i++)
	  {
		a=i;
		#pragma omp task firstprivate(i)//each task will have its own i determined when each task is created; each thread thus will have tasks with different i's.
		{
		  sleep(1);
		  printf("thread=%d, i=%d, a=%d\n", omp_get_thread_num(), i, a);//because a is shared, a=9 after all the tasks are created. at execution time, a=9.
		}
	  }
	}
	return 0;
}
