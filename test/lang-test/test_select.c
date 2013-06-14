#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* From Numerical Recipes*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
float Fselect_smallth(int k, int n, float arr[])
{/* to select the kth smallest value from an un-sorted array of n elements
  * WARNING: this function RE-ARRANGES arr to have arr[k-1] as the return value!
  * slightly modified from the original version to have arr[0~n-1] rather than arr[1~n]*/
	int i,ir,j,l,mid;
	float a,temp;

	arr--;//so that it's accessed with arr[1~n]
	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef SWAP

#define N 30
int main(int argc, char **argv)
{
	float arr[N]={21,5,2,15,10,29,6,27,16,13,7,4,28,20,24,30,26,25,18,14,22,1,3,17,23,12,8,19,9,11};
	float arr2[N];
	int k,i;
	k=atoi(argv[1]);
	memcpy(arr2,arr,N*sizeof(float));
	printf("%d,%d,%d\n",k,(int)Fselect_smallth(k,N,arr2),(int)arr2[k-1]);
	for(i=0;i<k-1;i++)
		printf("%d\t",(int)arr2[i]);
	printf("\n");
	
	for(i=0;i<k-1;i++)
		printf("%d\t",(int)arr[i]);
	printf("\n");
	
	return 0;
}
