#include <malloc.h>
#include <stdio.h>
#include <omp.h>

#define ORDER 1000
#define AVAL 3.0

/*
Compile the program:
$ gcc -fopenmp Recurso_M1_2b.c -o Recurso_M1_2b
*/

int main(int argc, char *argv[])
{
	int Ndim, Pdim;   /* A[N][P] */
	int i,j;
	
	double *A;
    double sum, serial_sum = 0.0;
	double start_time, run_time;


	Ndim = ORDER;
	Pdim = ORDER;

	A = (double *)malloc(Ndim*Pdim*sizeof(double));

	/* Initialize matrice */
	for (i=0; i<Ndim; i++)
		for (j=0; j<Pdim; j++) {
			*(A+(i*Ndim+j)) = AVAL;
			serial_sum += AVAL;
	}

	/* Do the sum of matrix elements */
	start_time = omp_get_wtime();
	//sum = 0.0; /* não é obrigatório igualar a 0*/
	#pragma omp parallel for reduction(+:sum) private(i, j)
	/*
	Alternativas:
	#pragma omp parallel for reduction(+:sum) collapse(2)
	*/
	for (i=0; i<Ndim; i++){
		for (j=0; j<Pdim; j++)
 			sum += *(A+(i*Ndim+j));
	}
	run_time = omp_get_wtime() - start_time;

	/* Check the answer */
	printf("Time: %f seconds \n", run_time);
    printf("%d threads\n",omp_get_max_threads());

    printf("serial_sum: %f\tparellel_sum: %f\n", serial_sum,sum);
	if (serial_sum != sum) 
		printf("Errors in sum\n");
	else
		printf("Hey, it worked\n");
}