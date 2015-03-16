#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int kk;

void paraprt(int n)
{
	//automatic vars in a routine called within an omp region is private implicitly
	int i,tid;
	//if you need shared vars inside omp function call, use static vars;
	//the side effect is static vars keeps its previous value if not reset a initial value during each call
	static int j,k;
	#pragma omp single
	{
	printf("tid: i,j\n");
	k=1;
	//~ j=10;
	}
	//~ j=k;
	printf("k=%d,j=%d\n",k,j);  //icc is dangerous!! gcc do it fine
	//~ printf("%d:k=%d,j=%d\n",omp_get_thread_num(),k,j);
	#pragma omp for reduction(+:k)
	for(i=0;i<n;i++)
	{
	tid=omp_get_thread_num();
	#pragma omp critical
	{
		j++;
		printf("%d: %d,%d\n",tid,i,j);
	}
	k++;
	}
	#pragma omp single
	kk=k;
}
int main()
{
	#pragma omp parallel num_threads(5)
	{
	paraprt(0);
	#pragma omp single
	printf("kk=%d\n",kk);
	paraprt(5);
	}
	printf("kk=%d\n",kk);

	return 0;
}
