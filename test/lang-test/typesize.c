#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//~ #define size_t long long
int main()
{
printf("short=%d\nint=%d\nlong int=%d\nlong long=%d\nfloat=%d\ndouble=%d\nlong double=%d\n",sizeof(short),sizeof(int),sizeof(long int),sizeof(long long),sizeof(float),sizeof(double),sizeof(long double));
printf("size_t=%d\n",sizeof(size_t));
printf("char=%d\n",sizeof(char));
return 0;
}
