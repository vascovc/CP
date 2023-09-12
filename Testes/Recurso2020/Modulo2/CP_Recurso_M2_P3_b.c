#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h> 

/*
Copile:
$ mpicc -o CP_Recurso_M2_P3_b CP_Recurso_M2_P3_b.c
Run:
$ mpirun -np 9 ./CP_Recurso_M2_P3_b
*/

int cmpfunc(const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main(int argc, char *argv[])
{	
	int myid, numprocs;
	int m, *rbuf, new_m;

	MPI_Init(&argc,&argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	m = (numprocs-1)-myid;
  	if (myid == 0)
		rbuf = (int *)malloc(numprocs*sizeof(int));

	MPI_Gather(&m, 1, MPI_INT, rbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (myid == 0) {
		for(int i = 0; i < numprocs; i++)
			printf("rbuf[%d]: %d\n",i,rbuf[i]);
		
		qsort(rbuf, numprocs, sizeof(int), cmpfunc);

		for(int i = 0; i < numprocs; i++)
			printf("rbuf_sorted[%d]: %d\n",i,rbuf[i]);
	}		
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* int MPI_Scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
					   void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) */
	MPI_Scatter(rbuf, 1, MPI_INT, &new_m, 1, MPI_INT, 0, MPI_COMM_WORLD);

	printf("P: %d | m: %d | new_m: %d\n", myid, m, new_m);

	MPI_Finalize();
  	return 0;
}