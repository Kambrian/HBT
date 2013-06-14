//feof do not tell the end when u just finished reading everything;
// it only tells the end when you have tried to move beyond EOF.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{
	FILE *fp;
	char a='a';
	fp=fopen("testeof.dat","w");
	fwrite(&a,sizeof(char),1,fp);
	fclose(fp);
	
    long n = 0, nread=0;
	fp=fopen("testeof.dat","r");
    while (!feof(fp)) 
	{
      nread+=fread (&a,sizeof(char),1,fp);
      n++;
    }
    fclose (fp);
    printf ("%d bytes, read %d times\n",nread,n);
  return 0;
}

