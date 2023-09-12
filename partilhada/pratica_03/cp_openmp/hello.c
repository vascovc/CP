// adapted from "A “Hands-on” Introduction to OpenMP" Tim Mattson Intel Corp. timothy.g.mattson@intel.com 

#include <stdio.h>
#include <omp.h>

int main ()  
{
    omp_set_num_threads(15); // alinea 2.4
    #pragma omp parallel //num_threads(20) // alinea 2.3
    {
    int ID = omp_get_thread_num();

    printf("Hello (%d) ", ID);
    printf("World (%d) \n", ID);
    }
}

// para o 2.2 e preciso criar um programa que chame o ./hello e fazendo o setenv