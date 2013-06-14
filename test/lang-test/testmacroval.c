/*if you haven't define a macro with a value, you cannot compare it with a integer, because it is just substitued to nothing and the comparison reads 
#if  ==1
which would produce an error
*/

#include <stdio.h>

//#define VAL
#define VAL 1

int main()
{
#if VAL==0
putchar('0');
#elif VAL==1
putchar('1');
#else
putchar('-');
#endif

return 0;
}
