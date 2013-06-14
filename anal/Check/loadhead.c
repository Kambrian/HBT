#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main()
{
	char snapdir[512]=SNAPSHOT_DIR;
	
	logfile=stdout;
	load_particle_header(99,snapdir);

	printf("%g,%g\n",headerA.mass[0],headerA.time);
	return 0;
}
