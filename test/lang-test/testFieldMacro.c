#include <stdio.h>
#include <stdlib.h>

#define NAME hey


int main()
{
struct
{
char hi[3];
char hey[4];
} say={"hi","hey"};

printf("%s\n",say.NAME);
return 0;
}
