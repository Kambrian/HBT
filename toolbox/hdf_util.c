#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdf_util.h"
// #include "mymath.h"

size_t load_hdfmatrixF(char *datafile, FloatMat var[], int nvar)
{//load a float32 matrix
hid_t    file_id;
herr_t      status;
size_t     i, j,nread=0;
H5T_class_t dataclass;
size_t element_size;

//~ printf("sizeof size_t, hsize_t: %u, %u\n",sizeof(size_t),sizeof(hsize_t));

/* open file from ex_lite1.c */
file_id = H5Fopen (datafile, H5F_ACC_RDONLY, H5P_DEFAULT);

for(i=0;i<nvar;i++)
{
	/* get the dimensions of the dataset */
status=H5LTget_dataset_ndims ( file_id,var[i].name, &var[i].dim);
if(status<0) 
{
  printf("Warning: dataset %s does not exist in %s\n", var[i].name, datafile);
  return 0;
}//failed to find the data
var[i].size=malloc(sizeof(hsize_t)*var[i].dim);
status = H5LTget_dataset_info(file_id,var[i].name,var[i].size,&dataclass,&element_size);

if(dataclass!=H5T_FLOAT||element_size!=sizeof(float))
{
	printf("error: dataset not float! %d, %zd\n", dataclass,element_size);
	exit(1);
}
/* allocate memory */
size_t nel;
nel=1;
//~ printf("dims: ");
for(j=0;j<var[i].dim;j++)
{
//~ printf("%lu * ",var[i].size[j]);
nel*=var[i].size[j];
}
//~ printf("\b\b  \n");
nread+=nel;
var[i].x=malloc(sizeof(float)*nel);

/* read dataset */
status = H5LTread_dataset_float(file_id,var[i].name,var[i].x);
//~ for(j=0;j<MIN(var[i].size[0],5);j++)
//~ printf("%g,",var[i].x[j]);
//~ puts("\n");
}

/* close file */
status = H5Fclose (file_id);

 return nread;
}

size_t load_hdfmatrixD(char *datafile, DoubleMat var[], int nvar)
{
hid_t    file_id;
herr_t      status;
size_t     i, j,nread=0;
H5T_class_t dataclass;
size_t element_size;

//~ printf("sizeof size_t, hsize_t: %u, %u\n",sizeof(size_t),sizeof(hsize_t));

/* open file from ex_lite1.c */
file_id = H5Fopen (datafile, H5F_ACC_RDONLY, H5P_DEFAULT);

for(i=0;i<nvar;i++)
{
	/* get the dimensions of the dataset */
status=H5LTget_dataset_ndims ( file_id,var[i].name, &var[i].dim);
if(status<0) 
{
  printf("Warning: dataset %s does not exist in %s\n", var[i].name, datafile);
  return 0;
}//failed to find the data
var[i].size=malloc(sizeof(hsize_t)*var[i].dim);
status = H5LTget_dataset_info(file_id,var[i].name,var[i].size,&dataclass,&element_size);

if(dataclass!=H5T_FLOAT||element_size!=sizeof(double))
{
	printf("error: dataset not double! %d, %zd\n", dataclass,element_size);
	exit(1);
}
/* allocate memory */
size_t nel;
nel=1;
//~ printf("dims: ");
for(j=0;j<var[i].dim;j++)
{
//~ printf("%lu * ",var[i].size[j]);
nel*=var[i].size[j];
}
//~ printf("\b\b  \n");
nread+=nel;
var[i].x=malloc(sizeof(double)*nel);

/* read dataset */
status = H5LTread_dataset_double(file_id,var[i].name,var[i].x);
//~ for(j=0;j<MIN(var[i].size[0],5);j++)
//~ printf("%g,",var[i].x[j]);
//~ puts("\n");
}

/* close file */
status = H5Fclose (file_id);

 return nread;
}

size_t load_hdfmatrix(char *datafile, GenericMat var[], int nvar, hid_t datatype)
{
	//datatype: H5T_NATIVE_FLOAT; H5T_NATIVE_DOUBLE;...
hid_t    file_id;
herr_t      status;
size_t     i, j,nread=0;

//~ printf("sizeof size_t, hsize_t: %u, %u\n",sizeof(size_t),sizeof(hsize_t));

/* open file from ex_lite1.c */
file_id = H5Fopen (datafile, H5F_ACC_RDONLY, H5P_DEFAULT);

for(i=0;i<nvar;i++)
{
	/* get the dimensions of the dataset */
status=H5LTget_dataset_ndims ( file_id,var[i].name, &var[i].dim);
if(status<0) 
{
  printf("Warning: dataset %s does not exist in %s\n", var[i].name, datafile);
  return 0;
}//failed to find the data
var[i].size=malloc(sizeof(hsize_t)*var[i].dim);
status = H5LTget_dataset_info(file_id,var[i].name,var[i].size,&var[i].dataclass,&var[i].element_size);

if(var[i].dataclass!=H5Tget_class(datatype)||var[i].element_size!=H5Tget_size(datatype))
{
	printf("error: datatype do not match, %d,%d; %zd,%zd\n",
				var[i].dataclass, H5Tget_class(datatype),
				var[i].element_size, H5Tget_size(datatype));
	exit(1);
}
/* allocate memory */
size_t nel;
nel=1;
//~ printf("dims: ");
for(j=0;j<var[i].dim;j++)
{
//~ printf("%lu * ",var[i].size[j]);
nel*=var[i].size[j];
}
//~ printf("\b\b  \n");
nread+=nel;
var[i].x=malloc(var[i].element_size*nel);

/* read dataset */
status = H5LTread_dataset(file_id,var[i].name,datatype,var[i].x);
}

/* close file */
status = H5Fclose (file_id);

 return nread;
}

void writeHDFmatrix(hid_t file, const void * buf, const char * name, hsize_t ndim, const hsize_t *dims, hid_t dtype, hid_t dtype_file)
{
  hid_t dataspace = H5Screate_simple (ndim, dims, NULL);
  hid_t dataset= H5Dcreate2(file, name, dtype_file, dataspace, H5P_DEFAULT, H5P_DEFAULT,
                H5P_DEFAULT);
  if(!(NULL==buf||0==dims[0]))
  {
	herr_t status = H5Dwrite (dataset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  }
  H5Sclose(dataspace);
  H5Dclose(dataset);
}

int GetDatasetDims(hid_t dset, hsize_t dims[])
{
  hid_t dspace=H5Dget_space(dset);
  int ndim=H5Sget_simple_extent_dims(dspace, dims, NULL);
  H5Sclose(dspace);
  return ndim;
}