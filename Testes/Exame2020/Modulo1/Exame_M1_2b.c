#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>
#include<time.h> 

/*
Compile the program:
$ gcc -fopenmp Exame_M1_2b.c -o Exame_M1_2b
*/

int main() 
{ 
	// This pointer will hold the 
	// base address of the block created 
	int *ptr; 
	int n, i, pos; 
	int lower = 0, upper = 10;

	// Get the number of elements for the array 
	n = 10; 
	printf("Number of elements: %d\n", n); 

	// Dynamically allocate memory using malloc() 
	ptr = (int*)malloc(n * sizeof(int)); 

	// Check if the memory has been successfully allocated by malloc or not 
	if (ptr == NULL) { 
		printf("Memory not allocated.\n"); 
		exit(0); 
	} 
	else {
        /* Allocating Random Number Values To The Elements Of An Array */
		srand(time(0)); // Use current time as seed for random generator
		for (i = 0; i < n; i++)
			ptr[i] = (rand() % (upper-lower+1))+lower;
  		
        // Print the elements of the array 
        printf("The elements of the array are: "); 
        for (i = 0; i < n; ++i) { 
            printf("%d ", ptr[i]); 
        }
        printf("\n"); 

        /* Serial Calculation */
		int max_val_serial = ptr[0];
		for (i = 1; i < n; i++) {
			if (ptr[i] > max_val_serial)
				max_val_serial = ptr[i];
		}

		/* Parallel Calculation */
        int max_val = ptr[0];
        #pragma omp parallel for reduction(max:max_val) firstprivate(n) //num_threads(4)
    	for (i = 1; i < n; i++) {
       		max_val = max_val > ptr[i] ? max_val : ptr[i];
       		/*
       		if(ptr[i] > max_val)
       			max_val = ptr[i];
       		*/
       		printf("%d) max_val= %d \t ptr[%d]=%d\n",omp_get_thread_num(),max_val,i,ptr[i]);
       	}

       	/* Checking For Output Validity */

		if (max_val_serial == max_val)
			printf("\nThe Max Value Is Same From Serial And Parallel OpenMP Directive\n");
		else
			printf("\nERROR: The Max Value Is Not Same In Serial And Parallel OpenMP Directive\n");
       	printf("max_val: %d\n",max_val); 
		// Free the memory 
		free(ptr);
	}
	return 0; 
} 