#include <stdlib.h>
#include <mcheck.h>

extern int* newmem();

int main() { 

int *b;
	mtrace(); /* Starts the recording of memory allocations and releases */
	b=newmem();
//~ free(b);
	muntrace();

	return 0; /* exit */

}
