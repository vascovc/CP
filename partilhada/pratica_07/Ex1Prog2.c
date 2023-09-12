#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 1000

int main(int argc, char **argv) {
    int rank, size, i, j, k;
    double start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize matrices
    int local_n = N / size;
    double *A = (double *)malloc(N * local_n * sizeof(double));
    double *B = (double *)malloc(N * N * sizeof(double));
    double *C = (double *)calloc(N * local_n, sizeof(double));

    // Fill matrices with random values
    for (i = 0; i < N * local_n; i++) {
        A[i] = (double)rand() / RAND_MAX;
    }
    for (i = 0; i < N * N; i++) {
        B[i] = (double)rand() / RAND_MAX;
    }

    // Start timer
    start_time = MPI_Wtime();

    // Perform matrix multiplication
    for (i = 0; i < local_n; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }

    // Gather results from all processes
    MPI_Allgather(C, N * local_n, MPI_DOUBLE, C, N * local_n, MPI_DOUBLE, MPI_COMM_WORLD);

    // End timer
    end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Time taken: %f seconds\n", end_time - start_time);
    }

    // Free memory
    free(A);
    free(B);
    free(C);

    MPI_Finalize();
    return 0;
}
