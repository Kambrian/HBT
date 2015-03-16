#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{
	FILE *fp;
	char a;
	fp=fopen("testwrite.dat","r");
	fread(&a,1,1,fp);
	fseek(fp,0,SEEK_CUR);
	printf("%ld\n",ftell(fp));
	fread(&a,1,1,fp);
	fseek(fp,1,SEEK_CUR);
	printf("%ld\n",ftell(fp));
	fclose(fp);
	return 0;
}
	
