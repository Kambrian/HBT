#include <stdio.h>
#include <stdlib.h>

void fun(int a)
{
	if(a<0)
	return;
	printf("%d\n",a);
}
int main()
{
fun(-1);
return 0;
}
