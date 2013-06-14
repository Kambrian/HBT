#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ printf("Error opening file '%s'\n",filename);	exit(1);	}

int main()
{
	FILE *fp;
	int tmp;
	long long int nid,i,j;
	//~ myfopen(fp,"longrecord.dat","r");
	if(!(fp=fopen("longrecord.dat","rb")))
	{
		printf("error opening longrecord.dat\n");
		exit(1);
	}
	printf("int=%d,long long=%d\n",sizeof(int),sizeof(long long int));
	nid=1024*1024*1024;
	i=0;
	do
	{
	fread(&tmp,sizeof(int),1,fp);
	i++;
	if(tmp)
	printf(" (%lld,%lld),%d,%#x,%f\n",i,i*4,tmp,tmp,*(float *)&tmp);
	}while(!feof(fp));
	printf("%lld-%lld-2=%lld\n",i,nid*3,i-nid*3-2);
	fclose(fp);
}
