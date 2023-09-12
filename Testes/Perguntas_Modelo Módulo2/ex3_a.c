#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/*
Copile:
$ mpicc -o ex3_a ex3_a.c
Run:
$ mpirun -np 8 ./ex3_a
*/

int main(int argc, char *argv[])
{
	int rank, size;
	int color;
	int subrank, subsize;
	MPI_Comm subcomm;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	color = rank%4;
	printf("rank: %d | color: %d\n", rank,color);

	MPI_Barrier(MPI_COMM_WORLD);
	usleep(1E6);
	if (rank == 0)
		printf("\n");

	/* int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) 
	
	This function partitions the group associated with comm into disjoint subgroups, one for
	each value of 'color'. 
	Each subgroup contains all processes of the same 'color'. 
	Within each subgroup, the processes are ranked in the order defined by the value of the 
	argument key, with ties broken according to their rank in the old group.*/
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &subcomm);
	MPI_Comm_rank(subcomm, &subrank); // o rank do processo no seu novo comuncador
	MPI_Comm_size(subcomm, &subsize); // o tamanho do novo comuncador

	printf("rank: %d | subrank: %d | subsize: %d\n", rank,subrank,subsize);

	MPI_Comm_free(&subcomm);
	MPI_Finalize();

	return 0;	
}