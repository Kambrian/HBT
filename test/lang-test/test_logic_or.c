#include <stdio.h>
#include <stdlib.h>

int cmp_xy(int x, int y)
{
if(x==y)
return 0;

if(x<0||x>y)
return 1;

if(y<0||x<y)
return -1;

return -1;
}

int main(int argc, char **argv)
{
int x,y;
x=atoi(argv[1]);
y=atoi(argv[2]);

printf("%d\n",cmp_xy(x,y));

return 0;
}
