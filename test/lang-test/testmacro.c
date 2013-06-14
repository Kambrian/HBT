#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define myprint fprintf

int main()
{
	myprint(stdout,"hello from my print %d\n",1);
	return 0;
}
