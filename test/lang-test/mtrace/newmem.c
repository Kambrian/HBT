#include <stdio.h>
#include <stdlib.h>

int* newmem()
{
	int* a = NULL;

	a = malloc(sizeof(int)); /* allocate memory and assign it to the pointer */
	return a;
	//~ free(a); /* we free the memory we allocated so we don't have leaks */
}
