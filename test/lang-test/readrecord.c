#include <stdlib.h>
#include <stdio.h>

#define NID 1073741824
#define NELE 3221225486

//~ int pid[NELE];
int pid2[3][NID];
extern size_t C_read_fortran_record(void *ptr,size_t bytes,FILE *fp);
int main()
{
	FILE *fp;
	long int i,j;
	
	if(!(fp=fopen("longrecord.dat","r")))
	{
	printf("error opening file\n");
	exit;
	}
	printf("sizeofint:%u,size_t:%d\n",sizeof(int),sizeof(size_t));

	C_read_fortran_record(&pid2[0][0],NID*sizeof(int)*3,fp);
	for(i=0;i<3;i++)
		for(j=0;j<NID;j++)
			if(pid2[i][j]!=j) {printf("%ld,%ld,%d,%0x\n",i,j,pid2[i][j],pid2[i][j]);break;}
	printf("file position:%ld\n",ftell(fp));
	fclose(fp);
	return 0;
}

#define SKIP fread(&dummy,4,1,fp)
#define SKIP2 fread(&dummy2,4,1,fp)
#define CHECK if(dummy!=dummy2){printf("error!record brackets not match:%d,%d\n",dummy,dummy2);fflush(stdout);exit(1);} 
//~ #define size_t long long
size_t C_read_fortran_record(void *ptr,size_t bytes,FILE *fp)
{
	short flag_subrec;
	int dummy,dummy2;
	char *p;
	size_t nread;
	
	p=(char *)ptr;
	if(sizeof(char)!=1){printf("wrong char size=%d\n",sizeof(char));exit(1);};
	flag_subrec=0; nread=0;
	SKIP;
	if(dummy<0)
	{
		dummy*=-1;
		flag_subrec=1;
	}
	fread(p,1,dummy,fp);
	nread+=dummy;
	SKIP2;
	if(dummy2<0)
	{
	printf("error reading fortran record! \nsubrecords before FIRST subrec:%d,%d\n",dummy,dummy2);
	return -1*dummy;
	}
	CHECK;
	p+=dummy;
	while(flag_subrec)
	{
		flag_subrec=0;
		SKIP;
		if(dummy<0) 
		{
		dummy*=-1;
		flag_subrec=1;
		}
		fread(p,1,dummy,fp);
		nread+=dummy;
		SKIP2;
		if(dummy2>0)
		{
		printf("error reading fortran record! \nFIRST subrec appear late?:%d,%d\n",dummy,dummy2);
		return -1*dummy;
		}
		dummy2*=-1;
		CHECK;
		p+=dummy;
	}
	if(nread!=bytes)
	{
	printf("reading failed: want %ld bytes, got %ld bytes\n",bytes,nread);
	return -1;
	}
	return nread;
}
