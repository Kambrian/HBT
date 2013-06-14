#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
char buf[1024],str[1024],sarr[1024];
strcpy(buf,"hello");
sprintf(str,"world");
sprintf(sarr,"%s","again");
printf("%s\n%s\n%s\n",buf,str,sarr);
return 0;
}
