// adapted from "A “Hands-on” Introduction to OpenMP" Tim Mattson Intel Corp. timothy.g.mattson@intel.com 

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifndef N
#define N 5
#endif
#ifndef FS
#define FS 40
#endif

struct node {
   int data;
   int fibdata;
   struct node* next;
};

int fib(int n) {
   int x, y;
   if (n < 2) {
      return (n);
   } else {
      x = fib(n - 1);
      y = fib(n - 2);
	  return (x + y);
   }
}

void processwork(struct node* p) 
{
   int n;
   n = p->data;
   p->fibdata = fib(n);
}

struct node* init_list(struct node* p) {
    int i;
    struct node* head = NULL;
    struct node* temp = NULL;
    
    head = malloc(sizeof(struct node));
    p = head;
    p->data = FS;
    p->fibdata = 0;
    for (i=0; i< N; i++) {
       temp  =  malloc(sizeof(struct node));
       p->next = temp;
       p = temp;
       p->data = FS + i + 1;
       p->fibdata = i+1;
    }
    p->next = NULL;
    return head;
}

int main(int argc, char *argv[]) {
     double start, end;
     struct node *p=NULL;
     struct node *p1=NULL;
     struct node *p2=NULL;
     struct node *temp=NULL;
     struct node *head=NULL;
     
	 printf("Process linked list\n");
     printf("  Each linked list node will be processed by function 'processwork()'\n");
     printf("  Each ll node will compute %d fibonacci numbers beginning with %d\n",N,FS);      
 
     p = init_list(p);
     head = p;

     start = omp_get_wtime();
     {
        while (p != NULL) {
		   processwork(p);
		   p = p->next;
        }
     }

     end = omp_get_wtime();
     p = head;
	 while (p != NULL) {
        printf("%d : %d\n",p->data, p->fibdata);
        temp = p->next;
        free (p);
        p = temp;
     }  
	 free (p);

     printf("Compute Time: %f seconds\n", end - start);

    printf("5.1\n");
    // 5.1
     printf("OpenMP without task.\n");
     p1= init_list(p1);
     head = p1;

     double start_51 = omp_get_wtime();
     {
         struct node *parr[N+1];
         parr[0] = p1;
         int i = 1;
        while (p1 != NULL) {
		   p1 = p1->next;
         parr[i]=p1;
         i++;
        }

         #pragma omp parallel for schedule(static,1)
         for (int i = 0; i < N+1; i++) {
            processwork(parr[i]);
         }
     }

     double end_51 = omp_get_wtime();
     p1 = head;
	 while (p1 != NULL) {
        printf("%d : %d\n",p1->data, p1->fibdata);
        temp = p1->next;
        free (p1);
        p1 = temp;
     }  
	 free (p1);

     printf("Compute Time: %f seconds\n", end_51 - start_51);
     printf("Speed up: %f seconds\n",(end-start)/(end_51 - start_51));

     // 5.2
     printf("5.2\n");
     p2 = init_list(p2);
     head = p2;

     double start_52 = omp_get_wtime();
     #pragma omp parallel
     {
        #pragma omp single
        {
        while (p2 != NULL) {
           #pragma omp task
           {
		    processwork(p2);
           }
		   p2 = p2->next;
        }
        }
     }

     double end_52 = omp_get_wtime();
     p2 = head;
	 while (p2 != NULL) {
        printf("%d : %d\n",p2->data, p2->fibdata);
        temp = p2->next;
        free (p2);
        p2 = temp;
     } 
	 free (p2);

     printf("Compute Time: %f seconds\n", end_52 - start_52);
     printf("Speedup: %f \n", (end - start)/(end_52 - start_52));
     return 0;
}
