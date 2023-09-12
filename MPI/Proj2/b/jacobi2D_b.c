// Projeto 2
//
// Alexandre Rodrigues 92993
// Gustavo Morais 92978

// b)

// TODO:
// Adapte o programa para usar condições fronteira periódicas em ambas as direções, 𝑥 e 𝑦. 
// Note que, com fronteiras periódicas, os valores 
// 𝑉𝑖,𝑗 nas fronteiras do domínio não são fixados à partida, e têm de ser calculados da mesma forma 
// que os do interior do domínio, recorrendo aos seus 4 vizinhos. Mais concretamente, os pontos de 
// índices  (0,𝑗)  e  (𝑛𝑦 −1,𝑗)  são  vizinhos,  assim  como  o  são  os  pontos  (𝑖,0)  e  (𝑖,𝑛𝑥 −1),  o  que 
// implicará comunicações adicionais  ‘através’  das fronteiras. Note que, devido às condições 
// fronteira periódicas, em cada linha e coluna o número de pontos é igual ao número de intervalos 
// entre  pontos,  o  que  significa  que  ℎ𝑥 =2/𝑛𝑥  e  ℎ𝑦 =2/𝑛𝑦.


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b))

#define NXMAX 500
#define TOL 1e-6
#define MAXIT 5e5
#define L 1

double f(double x, double y)
{
    return 7*sin(2*M_PI*x)*cos(3*M_PI*x)*sin(2*M_PI*y)*cos(3*M_PI*y);
}

int main(int argc, char *argv[])
{
    int nprocs;
    int myid; 
    int nx, ny;

    int manager_rank = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    if (myid == manager_rank)
    {
        printf("Introduza numero de pontos {max %d, 0 para sair}: \n",NXMAX);
        scanf(" %d", &nx);
    }
    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    ny = nx;

    if (nx == manager_rank)
    {
        MPI_Finalize();
        return 0;
    }

    if (nx < 0 || nx > NXMAX)
    {
        MPI_Finalize();
        return 1;
    }

    int nprocs_col = (int) nprocs/2;

    int ndims = 2;
    int dims[2] = {nprocs_col, 2};
    // Alterado para ativar as fronteiras periódicas
    int periodic[2] = {1,1};
    MPI_Comm comm2D;
    int newid;
    int nbrbottom, nbrtop, nbrleft, nbrright;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, 1, &comm2D);

    if (comm2D == MPI_COMM_NULL) // processo extra é terminado (se nprocs é impar)
    {
        MPI_Finalize();
        return 0;
    }

    nprocs = nprocs_col * dims[1];

    MPI_Comm_rank(comm2D, &newid);

    MPI_Cart_shift(comm2D, 0, 1, &nbrbottom , &nbrtop);
    MPI_Cart_shift(comm2D, 1, 1, &nbrleft , &nbrright);

    printf("myid=%d, newid=%d, bot=%d, top=%d, left=%d, right=%d\n", myid, newid, nbrbottom, nbrtop, nbrleft, nbrright);   

    int firstrow, firstcol;
    int myrows, mycols;

    if (newid == manager_rank)
    {   
        int listfirstrow[nprocs];
        int listmyrows[nprocs];

        int listfirstcol[nprocs];
        int listmycols[nprocs];

        int nrows = (int)((double)(ny-2)/(double)nprocs_col + 0.5);
        
        // Linhas
        // Alterado para incluir as linhas de fronteira
        for (int i = 0; i < nprocs_col; i++)
        {
            listfirstrow[2*i] = i *  (nrows);
            listmyrows[2*i] = nrows+1; 
            listfirstrow[2*i+1] = i *  (nrows);
            listmyrows[2*i+1] = nrows+1;
        }
        // Altera o numero de linhas do penultimo e do ultimo
        listmyrows[nprocs-2] = ny - (nprocs_col - 1) * nrows;
        listmyrows[nprocs-1] = ny - (nprocs_col - 1) * nrows;
        

        // Colunas
        // Alterado para incluir as colunas de fronteira
        int ncols_temp = (int)((nx-2)/2);
        for (int i = 0; i < nprocs_col; i++)
        {
            listfirstcol[2*i] = 0;
            listmycols[2*i] = ncols_temp + 1;
            listfirstcol[2*i+1] = ncols_temp + 1;
            listmycols[2*i+1] = nx - 1 - ncols_temp;
        }

        MPI_Scatter(listfirstrow, 1, MPI_INT, &firstrow, 1, MPI_INT, newid, comm2D);
        MPI_Scatter(listmyrows, 1, MPI_INT, &myrows, 1, MPI_INT, newid, comm2D);

        MPI_Scatter(listfirstcol, 1, MPI_INT, &firstcol, 1, MPI_INT, newid, comm2D);
        MPI_Scatter(listmycols, 1, MPI_INT, &mycols, 1, MPI_INT, newid, comm2D);

    }
    else
    {
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstrow, 1, MPI_INT, manager_rank, comm2D);
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &myrows, 1, MPI_INT, manager_rank, comm2D);

        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstcol, 1, MPI_INT, manager_rank, comm2D);
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &mycols, 1, MPI_INT, manager_rank, comm2D);
    }

    MPI_Barrier(comm2D);
    printf("newid=%d, firstrow=%d, lastrow=%d, firstcol=%d, lastcol=%d\n", newid, firstrow, firstrow+myrows-1, firstcol, firstcol+mycols-1);


    double (*Vold)[mycols+2], (*Vnew)[mycols+2], (*myf)[mycols+2];
    Vold = calloc(myrows + 2, sizeof(*Vold));
    Vnew = calloc(myrows + 2, sizeof(*Vnew));
    myf = calloc(myrows + 2, sizeof(*myf));

    double h = ((double)2 * L) / ((double) nx);

    for (int j = 1; j < mycols + 1 ; j++)
    {
        for (int i = 1; i < myrows + 1; i++)
        {
            myf[i][j] = f(-L + (firstcol + j - 1) * h, -L + (firstrow + i - 1) * h);
        }
        
    }

    // Os ciclos para definir os valores nas fronteiras foram removidos

    MPI_Datatype column;
    MPI_Type_vector(myrows + 2, 1, mycols + 2 , MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    
    double tm1 = MPI_Wtime();

    for (int iter = 0; iter < MAXIT; iter++)
    {
        double sums[2] = {0.0, 0.0};
        double global_sums[2];

        for (int j = 1; j < mycols + 1 ; j++)
        {
            for (int i = 1; i < myrows + 1; i++)
            {
                Vnew[i][j] = (Vold[i+1][j] + Vold[i-1][j] + Vold[i][j+1] + Vold[i][j-1]  - h * h  * myf[i][j]) / 4.0;
                sums[0] += (Vnew[i][j] - Vold[i][j]) * (Vnew[i][j] - Vold[i][j]);
                sums[1] += Vnew[i][j] * Vnew[i][j];
            }
            
        }

        MPI_Allreduce(sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, comm2D);
        
        if (sqrt(global_sums[0]/global_sums[1]) < TOL)
        {
            if (newid == manager_rank)
            {
                printf("Calculo demorou %f seg; %d iteracoes\n", MPI_Wtime()-tm1, iter);
                tm1 = MPI_Wtime();
            }

            // Agora nao é necessario ter em conta as colunas fantasmas
            // Isto torna o calculo de lsizes e start_ind mais fácil (sem verificações extras)
            int gsizes[2] = {ny, nx};
            int lsizes[2] = {myrows, mycols}; 
            int start_ind[2] = {firstrow, firstcol};

            MPI_Datatype filetype;
            MPI_Type_create_subarray(2, gsizes, lsizes, start_ind, MPI_ORDER_C, MPI_DOUBLE, &filetype);
            MPI_Type_commit(&filetype);
          
            int memsizes[2] = {myrows+2, mycols+2};
            start_ind[0] = 1;
            start_ind[1] = newid % 2;
            if (newid == 0 || newid == 1) {
                start_ind[0]--;
            }

            MPI_Datatype memtype;
            MPI_Type_create_subarray(2, memsizes, lsizes, start_ind, MPI_ORDER_C, MPI_DOUBLE, &memtype);
            MPI_Type_commit(&memtype);

            MPI_File fp;
            MPI_File_open(comm2D, "results_b.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
            MPI_File_set_view(fp, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
            
            MPI_File_write_all(fp, Vnew, 1, memtype, MPI_STATUS_IGNORE);
            MPI_File_close(&fp);

            MPI_Type_free(&filetype);
            MPI_Type_free(&memtype);

            if (newid == manager_rank)
            {
                printf("Escrita demorou %f seg\n", MPI_Wtime()-tm1);
            }

            break;
        }

        // comunicações sentido ascendente
        MPI_Sendrecv(Vnew[myrows], mycols+2, MPI_DOUBLE, nbrtop, 0, Vnew[0] , mycols+2, MPI_DOUBLE, nbrbottom, 0, comm2D, MPI_STATUS_IGNORE);
        
        // comunicações sentido descendente
        MPI_Sendrecv(Vnew[1], mycols+2, MPI_DOUBLE, nbrbottom, 1, Vnew[myrows+1] , mycols+2, MPI_DOUBLE, nbrtop, 1, comm2D, MPI_STATUS_IGNORE);
        
        // comunicações sentido para direita
        MPI_Sendrecv(&(Vnew[0][mycols]), 1, column, nbrright, 2, &(Vnew[0][0]), 1, column, nbrleft, 2, comm2D, MPI_STATUS_IGNORE);

        // comunicações sentido para esquerda
        MPI_Sendrecv(&(Vnew[0][1]), 1, column, nbrleft, 3, &(Vnew[0][mycols+1]), 1, column, nbrright, 3, comm2D, MPI_STATUS_IGNORE);
        
        for (int i = 0; i < myrows + 2; i++)
        {
            for (int j = 0; j < mycols + 2; j++)
            {
                Vold[i][j] = Vnew[i][j];
            }
            
        }
        
    }

    MPI_Type_free(&column);

    free(Vold);
    free(Vnew);
    free(myf);

    MPI_Finalize();

    return 0;
}