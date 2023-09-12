#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/*
Copile:
$ mpicc -o CP_Exame_M2_P4 CP_Exame_M2_P4.c
Run:
$ mpirun -np 3 ./CP_Exame_M2_P4
*/

int main (int argc, char *argv[]) {
	double a[10] = {1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0,};
	int b[10] = {1,2,3,4,5,6,7,8,9,10};
	int i, N;
	int myid, numprocs;
	int lsizes[1], gsizes[1], start_index[1];

	MPI_File fh;
  	MPI_Datatype filetype;

	MPI_Init(&argc,&argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
 	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

 	// 10 int ocupam tanto espaço como 5 doubles
	N = 10/2;

	gsizes[0] = 3*N;
	lsizes[0] = N;
	start_index[0] = myid*N;

	MPI_Type_create_subarray(1, gsizes, lsizes, start_index, MPI_ORDER_C, MPI_DOUBLE, &filetype);
	MPI_Type_commit(&filetype);

	MPI_File_open(MPI_COMM_WORLD, "output_v1.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	/* int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,
							 MPI_Datatype filetype, const char *datarep, MPI_Info info) */
	MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

	/* int MPI_File_write_all(MPI_File fh, const void *buf, int count,
							  MPI_Datatype datatype, MPI_Status *status) */
	if(myid == 0)
		MPI_File_write_all(fh, &a[0], N, MPI_DOUBLE, MPI_STATUS_IGNORE);
	if(myid == 1)
		MPI_File_write_all(fh, &a[N], N, MPI_DOUBLE, MPI_STATUS_IGNORE);
	if(myid == 2)
		MPI_File_write_all(fh, &b[0], 2*N, MPI_INT, MPI_STATUS_IGNORE);

	MPI_File_close(&fh);

	MPI_Finalize();
    exit(0);
}
