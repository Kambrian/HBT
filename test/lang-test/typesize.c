#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
//~ #define size_t long long
struct groupV4_header
{
  int Ngroups;
  int Nsubgroups;
  int Nids;
  int TotNgroups;
  int TotNsubgroups;
  int TotNids;
  int num_files;//long? no, but padding may exist.
  double time;
  double redshift;
  double HubbleParam;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  int flag_doubleprecision;//long? no, but padding may exist
};
int main()
{
printf("short=%zd\nint=%zd\nlong int=%zd\nlong long=%zd\nfloat=%zd\ndouble=%zd\nlong double=%zd\n",sizeof(short),sizeof(int),sizeof(long int),sizeof(long long),sizeof(float),sizeof(double),sizeof(long double));
printf("size_t=%zd\n",sizeof(size_t));
printf("char=%zd\n",sizeof(char));
printf("size_max=%zd\n", SIZE_MAX);
printf("array=%zd,%ld\n", (size_t)(sizeof(long)*3*12413832192LL), (long)(sizeof(long)*3*12413832192LL));

printf("groupheader=%zd\n", sizeof(struct groupV4_header));
return 0;
}
