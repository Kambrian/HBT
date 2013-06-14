//union or struct do not need to have a name when appearing inside another composite type, 
//and can be quoted directly with memeber name
#include <stdio.h>

typedef union
{
int shortid[2];
long long longid;
} HBTPID;

typedef union
{
	struct{
	int id0;
	int id1;
	};
	long long longid;
} HBTPID2;

int main()
{
HBTPID myid[2];
HBTPID2 myid2;
struct
{
	int i;
	HBTPID2;
} a;
myid[0].longid=0x0000000100000002;
myid[1].shortid[0]=1;
myid[1].shortid[1]=2;
myid2.id0=1;
myid2.id1=2;
printf("size: %u, %u\n",sizeof(myid[0].shortid[0]),sizeof(myid[0].longid));
printf("%d,%d,%#llx\n",myid[0].shortid[0],myid[0].shortid[1],myid[0].longid);
printf("%d,%d,%#llx\n",myid[1].shortid[0],myid[1].shortid[1],myid[1].longid);
printf("%d,%d,%#llx\n",myid2.id0,myid2.id1,myid2.longid);

a.i=1;
a.id0=1;
a.id1=2;
printf("%d,%d,%d,%#llx\n",a.i,a.id0,a.id1,a.longid);
return 0;

}
