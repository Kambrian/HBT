//when both positive and no overflow, unsigned and signed integers are exactly the same in memory
#include <stdio.h>
#include <stdlib.h>

int main()
{
unsigned int x;
int y;
x=123;
y=123;
FILE *fp;
printf("%X\n%X\n",x,y);
fp=fopen("unsigned_int.dat","w");
fwrite(&x,sizeof(unsigned int),1,fp);
fclose(fp);
fp=fopen("signed_int.dat","w");
fwrite(&y,sizeof(int),1,fp);
fclose(fp);
return 0;
}
