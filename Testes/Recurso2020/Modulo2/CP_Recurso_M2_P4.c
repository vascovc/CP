#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
/*
Copile:
$ mpicc -o CP_Recurso_M2_P4 CP_Recurso_M2_P4.c -lm
Run:
$ mpirun -np 9 ./CP_Recurso_M2_P4
*/

int main(int argc, char *argv[])
{
	int ndims = 2;
  	int periods[2] = {false, false};
  	//int periods[2] = {true, true};
  	int reorder = true;
  	int dims[2];
  	int myid, numprocs;
  	int newid, nbrbottom, nbrtop, nbrleft, nbrright, topright;
  	MPI_Comm comm2d;

	  MPI_Init(&argc,&argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myid);  	

  	/* Get a new communicator for a decomposition of the domain */
  	dims[0] = sqrt(numprocs);
  	dims[1] = dims[0];

  	/* int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[],
							const int periods[], int reorder, MPI_Comm *comm_cart) */
  	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm2d);
  	/* O ultimo vai embora se for impar. */
  	if (comm2d == MPI_COMM_NULL) {
    	MPI_Finalize();
    	exit(1);
  	}

  	/* Get my position in this communicator, and my neighbors */
  	MPI_Comm_rank(comm2d, &newid);
  	MPI_Cart_shift(comm2d, 0,  1, &nbrtop, &nbrbottom);
  	MPI_Cart_shift(comm2d, 1,  1, &nbrleft, &nbrright);

  	/* int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) 
	
	comm -> communicator with Cartesian structure (handle)
	rank -> rank of a process within group of comm (integer)
	maxdims -> length of vector coords in the calling program (integer)
	OUT coords -> integer array (of size ndims) containing the Cartesian coordinates of 
				  specified process (array of integers)
  	*/
  	int coords[2];
  	MPI_Cart_coords(comm2d, newid, 2, coords);
	
	/* int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]) 
	The function MPI_CART_GET() returns the Cartesian topology information that was associated
	with a communicator by MPI_CART_CREATE().

	IN comm -> communicator with Cartesian structure (handle)
	IN maxdims -> length of vectors dims, periods, and coords in the calling program (integer)
	
	OUT dims -> number of processes for each Cartesian dimension (array of integer)
	OUT periods -> periodicity (true/false) for each Cartesian dimension (array of logical)
	OUT coords -> coordinates of calling
	
	int outdims[2], outperiods[2], outcoords[2];
  	MPI_Cart_get(comm2d, 2, outdims, outperiods, outcoords);
	*/

  	/* int MPI_Sendrecv(const void *sendbuf, 
  						int sendcount, 
  						MPI_Datatype sendtype, 
  						int dest, 
  						int sendtag, 

  						void *recvbuf, 
  						int recvcount,
						MPI_Datatype recvtype, 
						int source, 
						int recvtag, 
						MPI_Comm comm,
						MPI_Status *status) */
  	
  	topright = MPI_PROC_NULL;
	MPI_Sendrecv(&nbrtop, 1, MPI_INT, nbrleft, 0,
				 &topright, 1, MPI_INT, nbrright, 0, comm2d, MPI_STATUS_IGNORE);
  	
  	printf("myid %d | newid %d | (%d, %d) | top %d | bottom %d | left %d | right %d | topright %d \n", 
  		myid, newid, coords[0], coords[1], nbrtop, nbrbottom, nbrleft, nbrright, topright);
  	
  	// printf("myid %d | newid %d | topright %d \n", myid, newid, topright);
  	fflush(stdout);
  	MPI_Barrier(comm2d);
  	usleep(5E5);

  	MPI_Finalize();
    exit(0);
}