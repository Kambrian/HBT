#ifndef HDF_HEADER_INCLUDED
	#ifdef HDF_V16
		/* hdf5 v1.6 */
		#include "H5LT.h"
	#else
		/* hdf5 v1.8 */
		#include "hdf5.h"
		#include "hdf5_hl.h"
	#endif


typedef struct 
{
	void *x; 
	int dim; //dimension
	hsize_t *size;  
	H5T_class_t dataclass;
	size_t element_size;
	char name[1024]; //name of matrix
} GenericMat;

typedef struct 
{
	float *x; 
	int dim; //dimension
	hsize_t *size;  
	char name[1024]; //name of matrix
} FloatMat;

typedef struct 
{
	double *x; 
	int dim; //dimension
	hsize_t *size;  
	char name[1024]; //name of matrix
} DoubleMat;

#define FREE_HDFMATRIX(M) {free(M.size);free(M.x);}

extern size_t load_hdfmatrixF(char *datafile, FloatMat var[], int nvar);
extern size_t load_hdfmatrixD(char *datafile, DoubleMat var[], int nvar);
extern size_t load_hdfmatrix(char *datafile, GenericMat var[], int nvar, hid_t datatype);
//datatype: H5T_NATIVE_FLOAT; H5T_NATIVE_DOUBLE;...

#ifdef HDF_V16
#define HDFcreate_group(file,group) H5Gcreate1(file,group,16)
#else
#define HDFcreate_group(file,group) H5Gcreate2(file,group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

#define HDF_HEADER_INCLUDED
#endif
