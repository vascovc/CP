#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int num_procs, rank,h,n;

    double pi = 0.0;
    //double dx = 0.00001; // Step size
    //int n = 1000000000; // Number of intervals

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    while(1){

        if (rank==0){
            printf("Introduza o numero de intervalos: ");
            scanf(" %d",&n);
        }
        
        MPI_Bcast(&n,1,MPI_INTEGER,0,MPI_COMM_WORLD);

        if (n==0){
            break;
        }
        
        h = 1.0 /(double) n;
        double local_pi = 0.0;
        //for (int i = rank * local_n; i < (rank + 1) * local_n; i++) {
        for (int i=rank;i<n;i+=num_procs){
            //double x = (i + 0.5) * dx;
            double x = h*(i+0.5);
            local_pi += 4.0 / (1.0 + x * x);
        }

        MPI_Reduce(&local_pi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            //pi *= dx;
            pi *= h;
            pi /= 2; 
            printf("pi obtido  = %f\n", pi);
            printf("pi teorico = %f\n", M_PI);
        }
    }

    MPI_Finalize();

    return 0;
}
