#include <stdio.h>
#include <stdlib.h>

extern void open_fortran_file_(char *filename,int *fileno,int *endian,int *error);
extern void read_fortran_record1_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record2_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record4_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record8_(void *arr,long int *arr_len,int *fileno);
extern void close_fortran_file_(int *fileno);

int main()
{
	char filename[1024]="testwrite10.fdat";
	int i,fileno=10;
	int arr[30]={1},endian=0,error_no;
	long int arrsize=30;
	float f[2];
	printf("short %d\n",sizeof(short));
	open_fortran_file_(filename,&fileno,&endian,&error_no);
	read_fortran_record4_(arr,&arrsize,&fileno);
	arrsize=2;
	read_fortran_record4_(f,&arrsize,&fileno);
	for(i=0;i<30;i++)
	printf("%d\t",arr[i]);
	printf("\n");
	printf("%f,%f\n",f[0],f[1]);
	close_fortran_file_(&fileno);
	return 0;
}
