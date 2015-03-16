#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
//~ #include <omp.h>

int main(int argc,char **argv)
{
	int i;
	double a,b,c,d,e,f;

	a=b=c=d=e=f=0.;	
	#pragma omp parallel for reduction(+: a,b,c,d,e,f)
	for(i=0;i<10;i++)
	{
	a+=1.;
	b+=1.;
	c+=1.;
	d+=1.;
	e+=1.;
	f+=1.;
    //~ printf("i=%d:%d,  %g,%g,%g,%g,%g,%g\n",i,omp_get_thread_num(),a,b,c,d,e,f);
	}
	printf("%g,%g,%g,%g,%g,%g\n",a,b,c,d,e,f);
	return 0;
}
