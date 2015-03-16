#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>

int main()
{
	int i;
	double s;
	char *outs;
	time_t time0,time1;
	clock_t ctime0,ctime1,ct0,ct1;
	struct tms t0,t1;
	time0=time(NULL);
	outs=ctime(&time0);
	ctime0=clock();
	ct0=times(&t0);
	s=1.;
	while(clock()<ctime0+30*CLOCKS_PER_SEC)
	{
		s++;
	}
	printf("%g\n",s);
	clock_t ticks_per_sec=sysconf(_SC_CLK_TCK);
	time1=time(NULL);
	ctime1=clock();
	ct1=times(&t1);
	printf("returning\n%s",outs);
	printf("walltime:%ld sec \n", (time1-time0));
	printf("%ld clocks per sec\n", CLOCKS_PER_SEC);
	//~ #define CLOCKS_PER_SEC 1
	printf("cputime: %ld sec (%ld - %ld)\n", (ctime1-ctime0)/CLOCKS_PER_SEC, ctime0, ctime1);
	printf("cputimes: %ld sec (%ld - %ld)\n", (ct1-ct0)/ticks_per_sec, ct0,ct1);
	printf(" usertimes: %ld sec (%ld - %ld)\n", (t1.tms_utime-t0.tms_utime)/ticks_per_sec, t0.tms_utime, t1.tms_utime);
	printf("systimes: %ld sec (%ld - %ld )\n", (t1.tms_stime-t0.tms_stime)/ticks_per_sec, t0.tms_stime,t1.tms_stime);
	printf("child usertimes: %ld sec (%ld - %ld )\n", (t1.tms_cutime-t0.tms_cutime)/ticks_per_sec, t0.tms_cutime,t1.tms_cutime);
	printf("child systimes: %ld sec (%ld - %ld )\n", (t1.tms_cstime-t0.tms_cstime)/ticks_per_sec, t0.tms_cstime,t1.tms_cstime);
	return 0;
}
