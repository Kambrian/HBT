#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define MSnap MaxSnap

int main()
{
	int i;
	float a[MSnap],z[MSnap];
	char buf[1024];
	FILE *fp;
	
	logfile=stdout;
	
	for(i=IniSnap;i<MSnap;i++)
	{
		load_particle_header(i,SNAPSHOT_DIR);
		a[i]=header.time;
		//~ z[i]=header.redshift;
	}
		
	sprintf(buf,"%s/Redshift.dat",SUBCAT_DIR);
	myfopen(fp,buf,"w");
	for(i=IniSnap;i<MSnap;i++)
		fprintf(fp,"%d\t%g\n",i,a[i]);
	
	fclose(fp);	
	return 0;
}
