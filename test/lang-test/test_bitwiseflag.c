//encode flags bitwise
#include <stdio.h>
#include <stdlib.h>

#define get_bit(x,k) (((x)&(1<<k))>>k)
int main()
{
unsigned char x=0b00000110,f1,f2,f3;
f1=get_bit(x,0);
f2=get_bit(x,1);
f3=get_bit(x,2);
printf("%d,%d,%d\n",f1,f2,f3);
printf("%d,%d,%d\n",x&0x1,(x&0b10)>>1,(x&0b100)>>2);
return 0;
}
