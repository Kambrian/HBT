//GOD: format string matters!!
//Integer Promotion: for variable argument functions, all the arguments are pushed into a stack, 
//converted to some basic datatypes before pushing.
//on 32bit system, the basic integer typer is int32, so any short or char are converted to int32 first 
//and then stored in the stack, for later print. larger integers, like int64, have to take up 64bits
//the compiler find the datatype according to format string, then find the extension of the data in the stack
//according to the datatype.
//on 64bit systems, probably the base int is int64, so you do not have a problem, luckily.
#include <stdio.h>

int main()
{
int i=10;
long j=100;
unsigned k=3;
long long l=150;

printf("%lld\n",i);
printf("%lld\n",(long long)i);
printf("%ld,%lld,%lld\n",i,i,(long long)i);
printf("%ld,%lld,%lld\n",(long)i,(long long)i,(long long)i);
//~ printf("%ld,%ld\n",j,j);
//~ printf("%ld,%ld\n",k,k);
printf("%03d,%03d,%03d\n",(int)l,(int)l,i);
//~ printf("%u,%lu\n",k,k);
printf("%d\n",(int)l);

for(l=1;l<10240;l*=10)
{
	printf("LINE%03d/FILE.%03d\n",l,l);
}

}
