#include <stdio.h>
#include <stdlib.h>

int main()
{
int *p,*pp;

p=malloc(0);
printf("malloc(0)==>non-null: %x\n",p);

p=realloc(p,0);
printf("realloc(malloc(0),0)==>0: %x\n",p);

pp=NULL;
pp=realloc(pp,0);
printf("realloc(null,0)==>malloc(0)==>non-null: %x\n",pp);

pp=malloc(1);
pp=realloc(pp,0);
printf("realloc(malloc(1),0)==>0: %x\n",pp);

free(p);
free(pp);

return 0;
}
