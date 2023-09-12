#include <stdio.h>
#include <time.h>

#define SIZE   10 //(1000*4)
#define REPEAT (100000)

/**
 * sumarray using mmx instructions 
 * */
void sumarray_mmx( int *a, int *b, int *c, int size )
{

  for (int i=0;i<size;i+=2) {
    if (i + 2 <= size)
    {  __asm__ volatile
        ( // instruction         comment          
        "\n\t movq     %1,%%mm0     \t#"            // move a[i] para mm0
        "\n\t movq     %2,%%mm1     \t#"            // move b[i] para mm1
        //"\n\t paddusb    %%mm0,%%mm1    \t#"
        "\n\t paddd    %%mm0,%%mm1    \t#"          // faz mm1 = mm0 + mm1 em 2 blocos de 32 bits
        "\n\t movq     %%mm1,%0     \t#"            // move mm1 para c[i]
        : "=m" (c[i])      // %0
        : "m"  (a[i]),     // %1 
          "m"  (b[i])      // %2
        );  
    }
    else c[size-1] = a[size-1] + b[size-1];
  }

   __asm__("emms" : : );
}
/**
 * summarray using sse instructions
 * */

void sumarray_sse(int *a, int *b, int *c, int size)
{
    for (int i=0;i<size;i+=4) 
    {
        if (i + 4 <= size)
        {
            __asm__ volatile
                (
                "\n\t movdqa     %1,%%xmm0     \t#"
                "\n\t movdqa     %2,%%xmm1     \t#"
                "\n\t paddd    %%xmm0,%%xmm1    \t#" // 4 double words de 32 bits por registo
                "\n\t movdqa     %%xmm1,%0     \t#"
                : "=m" (c[i])      // %0
                : "m"  (a[i]),     // %1
                "m"  (b[i])      // %2
                );
        }
        else
        {
            for (int j = i; j < size; j++)
            {
                c[j] = a[j] + b[j];
            }
        }
    }
}

/**
 * sumarray using classic code 
 * */
void sumarray( int *a, int *b, int *c, int size )
{
  for (int i=0;i<size;i++) {
      c[i]=a[i]+b[i];
  }
}

/**
 * print array
 * */
void print_array(int *a, int size)
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
void initArrays( int *a, int *b, int *c, int size )
{
    for (int i=0; i< SIZE; i++) {
        a[i]=(i<<16)+1;
        b[i]=0xffff;
        c[i]=0;
    }
}


/**
 * test summation functions
 */
int main(void)
{
    int a[SIZE];
    int b[SIZE],c[SIZE];

    int n, nelemsum;

    clock_t init, end;

    //initialize arrays
    nelemsum=SIZE;
    initArrays(a,b,c,nelemsum);

    // test classic code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray(a,b,c,nelemsum);
    end = clock();

    print_array(c,SIZE+1);

    printf("sumarray time = %f\n", (end-init)/(CLOCKS_PER_SEC*1.0));

    //initialize arrays
    initArrays(a,b,c,nelemsum);

    // test mmx code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_mmx(a,b,c,nelemsum);
    end = clock();

    print_array(c,SIZE+1);

    printf("sumarray time = %f\n", (end-init)/(CLOCKS_PER_SEC*1.0));


    //initialize arrays
    initArrays(a,b,c,nelemsum);

    // test mmx code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_sse(a,b,c,nelemsum);
    end = clock();

    print_array(c,SIZE+1);

    printf("sumarray time = %f\n", (end-init)/(CLOCKS_PER_SEC*1.0));

    printf("\n");

    return 0;
}