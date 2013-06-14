#include <stdio.h>
#include <mpi.h>

MPI_Comm HBT_COMM_SLAVE;
int TaskID,SlaveID,NTask,NSlave;

void HBT_init_mpi(int argc, char **argv)
{
int master_range[3]={0,NMASTER-1,1};	
MPI_Group worldgroup, slavegroup;	
MPI_Comm MPI_COMM_SLAVE;

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&TaskID);
MPI_Comm_size(MPI_COMM_WORLD,&NTask);

MPI_Comm_group(MPI_COMM_WORLD,&worldgroup);
MPI_Group_range_excl(worldgroup,1,&master_range,slavegroup);
MPI_Comm_create(MPI_COMM_WORLD,slavegroup,&HBT_COMM_SLAVE);

MPI_Comm_rank(HBT_COMM_SLAVE,&SlaveID);
MPI_Comm_size(HBT_COMM_SLAVE,&NSlave);

MPI_Group_free(&worldgroup);
MPI_Group_free(&slavegroup);

printf("%d/%d, %d/%d\n",TaskID,NTask,SlaveID,NSlave);
}

int main(int argc, char ** argv)
{
HBT_init_mpi();
printf("\t %d/%d, %d/%d\n",TaskID,NTask,SlaveID,NSlave);
MPI_Finalize();
return 0;
}
