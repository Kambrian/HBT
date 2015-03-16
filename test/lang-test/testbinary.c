#include <stdio.h>
#include <stdlib.h>
/*==the binary files can be ported from this machine to altix===*/
int main()
{	
	FILE *fp;
	float pi=3.14159;

	fp=fopen("testwrite.dat","w");
	fwrite(&pi,sizeof(float),1,fp);
	fclose(fp);
	printf("%u\n",sizeof(float));
	fp=fopen("testwrite.dat","r");
	fread(&pi,sizeof(float),1,fp);
	fclose(fp);
	printf("%g\n",pi);
	return 0;
}
