#include <stdio.h>
#include <stdlib.h>
struct com
	{
		int a;
		int b;
	} ;
struct com * fun(struct com *x)
{
	return x+1;
}

int main()
{
	struct com x[2],*y;
	x[0].a=1;
	x[0].b=2;
	y=fun(x);
	y->a=3;
	fun(x)->b=4;
	printf("%d\n",x[0].a+fun(x)->a);
	return 0;
}
