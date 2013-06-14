#include <stdio.h>
#include <stdlib.h>

int main()
{
	int i,j;
	j=10;
	for(i=0;i<j;i++)
	{
		printf("%d,%d,%d,%d\n",i,j,--j,j);//decrease j before using it
	}
}
