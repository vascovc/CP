#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define NROW 100
#define NCOLS 100
#define min(a,b) (a<b?a:b)

int main(int argc, char *argv[]){

	int manager_rank=0;
	int b[NCOLS];
	int nprocs, my_rank;

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Status status;
	int tag_term = NROW;
	int c[NROW];
	if(my_rank == manager_rank){
		int A[NROW][NCOLS];
		// A = [0 1 2 3....; ncols ncols+1 ncols+2 ...]
		for(int i=0;i<NROW;i++){
			for(int j=0;j<NCOLS;j++){
				A[i][j] = j+i*NCOLS;
			}
		}

		for(int i=0;i<NCOLS;i++){
			b[i] = i;
		}

		MPI_Bcast(b,NCOLS,MPI_INT,manager_rank,MPI_COMM_WORLD); // broadcast ao vetor b
		int next_row=0;
		for(int i=1;i<=min(nprocs-1,NROW);i++){
			MPI_Send(A[next_row],NCOLS,MPI_INT,i,next_row,MPI_COMM_WORLD);
			next_row++;
		}

		int ans;
		for (int i = 0; i < NROW; i++)
		{
			MPI_Recv(&ans,1,MPI_INT,MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			int sender = status.MPI_SOURCE;
			int row_index = status.MPI_TAG;

			c[row_index] = ans;

			if(next_row < NROW) {
				MPI_Send(A[next_row],NCOLS,MPI_INT,sender,next_row,MPI_COMM_WORLD);
				next_row++;
			}else{
				MPI_Send(MPI_BOTTOM,0,MPI_INT,sender,tag_term,MPI_COMM_WORLD);
			}
		}
		for(int i =0; i<NROW; i++){
            printf("C[%d] = %d \n", i,c[i]);
        }
	}else{
		int row[NCOLS];
		MPI_Bcast(b,NCOLS,MPI_INT,manager_rank,MPI_COMM_WORLD);

		if(my_rank <= NROW){
			while(1){
				MPI_Recv(row,NCOLS,MPI_INT,manager_rank,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				int row_index = status.MPI_TAG;
				
				if(status.MPI_TAG == tag_term)
					break;

				int ans = 0;
				for(int i=0;i<NCOLS;i++){
					ans += row[i]*b[i];
				}
				MPI_Send(&ans,1,MPI_INT,manager_rank,row_index,MPI_COMM_WORLD);
		}
		}
	}
	printf("Process %d finished\n",my_rank);
	MPI_Finalize();
	
	return 0;
}
