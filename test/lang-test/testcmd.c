#include <stdio.h>
#include <stdlib.h>

int main()
{
int i=system("wc -l teststr.c");
printf("%d\n",i);
return 0;
}
