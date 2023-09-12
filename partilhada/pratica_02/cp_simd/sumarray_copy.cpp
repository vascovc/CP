
#include <stdio.h>
#include <time.h>

#define SIZE   (1000*4)
#define REPEAT (1000000)

/**
 * sumarray using mmx instructions 
 * */
void sumarray_mmx( char *a, char *b, char *c, int size )
{

  for (int i=0;i<size;i+=8) {
    __asm__ volatile
        ( // instruction         comment          
        "\n\t movq     %1,%%mm0     \t#"
        "\n\t movq     %2,%%mm1     \t#"
        "\n\t paddd    %%mm0,%%mm1    \t#" // o normal que ja vinha
        //"\n\t paddusb    %%mm0,%%mm1    \t#" // alinea 1.3 unsigned saturation byte //origina valores a faltar uma unidade devido a dar 
        "\n\t movq     %%mm1,%0     \t#"
        : "=m" (c[i])      // %0
        : "m"  (a[i]),     // %1 
          "m"  (b[i])      // %2
        );  
  }

   __asm__("emms" : : );
}

void sumarray_sse( char *a, char *b, char *c, int size )
{

  for (int i=0;i<size;i+=16) {
    __asm__ volatile
        ( // instruction         comment          
        "\n\t movdqa    %1,%%xmm0     \t#"
        "\n\t movdqa    %2,%%xmm1     \t#"
        "\n\t paddd    %%xmm0,%%xmm1    \t#" // o normal que ja vinha
        //"\n\t paddusb    %%mm0,%%mm1    \t#" // alinea 1.3 unsigned saturation byte //origina valores a faltar uma unidade devido a dar 
        "\n\t movdqa     %%xmm1,%0     \t#"
        : "=m" (c[i])      // %0
        : "m"  (a[i]),     // %1 
          "m"  (b[i])      // %2
        );  
  }
}

/**
 * sumarray using classic code 
 * */
void sumarray( char *a, char *b, char *c, int size )
{
  for (int i=0;i<size;i++) {
      c[i]=a[i]+b[i];
  }
}

/**
 * print array
 * */
void print_array(char *a, int size)
{
    printf("base10: ");
    for (int i=0; i < size; i++) {
    printf("%10d",a[i]);
    }
    printf("\nbase16: ");
    for (int i=0; i < size; i++) {
    printf("%10x",a[i]);
    }
    printf("\n");
}

/**
 * init arrays
 * */
void initArrays( char *a, char *b, char *c, int size )
{
    for (int i=0; i< SIZE; i++) {
        a[i]=(i<<4)+1;
        b[i]=0xff;
        c[i]=0;
    }
}


/**
 * test summation functions
 */
int main(void)
{
    //int a[SIZE];
    //int b[SIZE];
    //int c[SIZE];
    //para a 1.4
    char a[SIZE] __attribute__((aligned (16)));
    char b[SIZE] __attribute__((aligned (16)));
    char c[SIZE] __attribute__((aligned (16)));

    int n, nelemsum;

    clock_t init_seq, end_seq, init_mmx, end_mmx, init_sse, end_sse;
    
    //initialize arrays
    nelemsum=SIZE;
    initArrays(a,b,c,nelemsum);
    
    // test classic code
    init_seq = clock();
    for(n=0;n<REPEAT;n++)
        sumarray(a,b,c,nelemsum);
    end_seq = clock();

    print_array(c,12);
    
    float t_diff_seq = end_seq-init_seq;
    printf("sumarray time = %f\n", (t_diff_seq)/(CLOCKS_PER_SEC*1.0));

    //initialize arrays
    initArrays(a,b,c,nelemsum);
    
    // test mmx code
    init_mmx = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_mmx(a,b,c,nelemsum);
    end_mmx = clock();

    print_array(c,12);
    
    float t_diff_mmx = end_mmx-init_mmx;
    printf("sumarray time = %f\n", (t_diff_mmx)/(CLOCKS_PER_SEC*1.0));
    
    // speedup
    printf("speedup: %f\n", t_diff_seq/t_diff_mmx);

    //initialize arrays
    initArrays(a,b,c,nelemsum);
    
    // test sse code
    init_sse = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_sse(a,b,c,nelemsum);
    end_sse = clock();

    print_array(c,12);
    
    float t_diff_sse = end_sse-init_sse;
    printf("sumarray time = %f\n", (t_diff_sse)/(CLOCKS_PER_SEC*1.0));
    
    // speedup
    printf("speedup: %f\n", t_diff_mmx/t_diff_sse);

    printf("\n");

    return 0;
}

    
