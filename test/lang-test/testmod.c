#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
int i,j;
i=atoi(argv[1]);
j=atoi(argv[2]);
printf("%d\n",i%j);
return 0;
}
