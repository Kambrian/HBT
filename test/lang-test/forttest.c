#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

extern void write_fortran_file_(char *filename, int *fileno);
extern void read_fortran_file_(char *filename, int *fileno, int *i);
extern void open_fortran_file_(char *filename,int *fileno,int *endian,int *error);
extern void read_fortran_record1_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record2_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record4_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record8_(void *arr,long int *arr_len,int *fileno);
extern void close_fortran_file_(int *fileno);

int main()
{
	int i,fileno=11;
	
#pragma omp parallel//for num_threads(3)
#pragma omp single
	for(i=0;i<10;i++)
#pragma omp task firstprivate(i) shared(fileno)
	{
	  char buf[1024];
	  sprintf(buf, "file%d.dat", i);
	  int j=1,cfile=fileno+i, endian=0, error;
	  long int len=1;
// 	  printf("%d-%d:%d\n", i,omp_get_thread_num(), cfile);
// 	  write_fortran_file_(buf, &cfile);
// #pragma omp critical
	  open_fortran_file_(buf, &cfile, &endian,  &error);
// 	  read_fortran_record4_(&j, &len, &cfile);
	  close_fortran_file_(&cfile);
// #pragma omp critical
// 	  read_fortran_file_(buf, &cfile, &j);
// #pragma omp critical
	  printf("%d-%d:%d, %d\n",i, omp_get_thread_num(), cfile, j);
	}
	return 0;
}
