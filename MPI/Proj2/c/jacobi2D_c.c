// Projeto 2
//
// Alexandre Rodrigues 92993
// Gustavo Morais 92978

// c)

// TODO:
// A aproximação por diferenças finitas usada na alínea anterior para as segundas derivadas 
// corresponde a um estêncil de 5 pontos (o ponto (𝑖,𝑗) mais 4 vizinhos), que introduz um erro local 
// da ordem de ℎ2. A ordem deste erro pode ser reduzida para ℎ4 usando um estêncil de 9 pontos, 
// baseado em 𝜕2𝑉(𝑥,𝑦)
// 𝜕𝑥2 ≅−𝑉(𝑥−2ℎ,𝑦)+16 𝑉(𝑥−ℎ,𝑦)−30𝑉(𝑥,𝑦)+16 𝑉(𝑥+ℎ,𝑦)−𝑉(𝑥+2ℎ,𝑦)
// 12 ℎ2  para a derivada em 𝑥, 
// e na aproximação correspondente para a derivada em 𝑦. Este estêncil de  9 pontos usa  4 pontos 
// adicionais  a  uma  distância  de  2ℎ.  Quando  aplicado  ao  método  de  Jacobi  ponderado  (aqui,  a 
// ponderação traz estabilidade ao método), este estêncil produz a seguinte equação iterativa: 
// 𝑉𝑖,𝑗(𝑘) = 𝑤
// 60[16𝑉𝑖−1,𝑗(𝑘−1) +16𝑉𝑖+1,𝑗(𝑘−1) +16𝑉𝑖,𝑗−1(𝑘−1) +16𝑉𝑖,𝑗+1(𝑘−1) −𝑉𝑖−2,𝑗(𝑘−1)
// −𝑉𝑖+2,𝑗(𝑘−1) −𝑉𝑖,𝑗−2(𝑘−1) −𝑉𝑖,𝑗+2(𝑘−1) −12ℎ2𝑓𝑖,𝑗]+(1−𝑤)𝑉𝑖,𝑗(𝑘−1) . 
// Modifique  o  programa  desenvolvido  na  alínea  b)  de  modo  a  aplicar  este  algoritmo,  com  o 
// parâmetro  de  ponderação  𝑤 =15/16.  Repare  que,  agora,  cada  processo  tem  de  receber  duas 
// linhas  e  duas  colunas  dos  seus  processos  vizinhos.  Faça  a  representação  gráfica  do  resultado  e 
// compare-a com a solução obtida na alínea anterior. 

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b))

#define NXMAX 500
#define TOL 1e-6
#define MAXIT 5e5
#define L 1.0

#define W (15.0/16.0)

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

    MPI_Bcast(&nx, 1, MPI_INT, manager_rank, MPI_COMM_WORLD);
    ny = nx;

    if (nx == 0)
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
    
    MPI_Barrier(comm2D);

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

    // ALterado para usar mais 2 colunas e linhas fantasma, necesário para o maior estêncil
    double (*Vold)[mycols+4], (*Vnew)[mycols+4], (*myf)[mycols+4];
    Vold = calloc(myrows + 4, sizeof(*Vold));
    Vnew = calloc(myrows + 4, sizeof(*Vnew));
    myf = calloc(myrows + 4, sizeof(*myf));

    double h = ((double)2 * L) / ((double) nx);

    // dominio principal agora nao inclui as 2 primeiras e 2 ultimas linhas e colunas
    for (int j = 2; j < mycols + 2 ; j++)
    {
        for (int i = 2; i < myrows + 2; i++)
        {
            myf[i][j] = f(-L + (firstcol + j - 2) * h, -L + (firstrow + i - 2) * h);
        }    
    }

    // Alterado para + 4
    MPI_Datatype column;
    MPI_Type_vector(myrows + 4, 1, mycols + 4, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    
    double tm1 = MPI_Wtime();

    for (int iter = 0; iter < MAXIT; iter++)
    {
        double sums[2] = {0.0, 0.0};
        double global_sums[2];

        // Dominio principal agora nao inclui as 2 primeiras e 2 ultimas linhas e colunas
        // Alterado para usar o novo estêncil e nova equação iterativa
        for (int j = 2; j < mycols + 2 ; j++)
        {
            for (int i = 2; i < myrows + 2; i++)
            {
                Vnew[i][j] = (W/60)*(16*Vold[i-1][j] + 16*Vold[i+1][j] + 16*Vold[i][j-1]+ 16*Vold[i][j+1] 
                     - Vold[i-2][j] - Vold[i+2][j] - Vold[i][j-2] - Vold[i][j+2] -12*h*h*myf[i][j]) + (1-W)*Vold[i][j];
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

            int gsizes[2] = {ny, nx};
            int lsizes[2] = {myrows, mycols};
            int start_ind[2] = {firstrow, firstcol};

            MPI_Datatype filetype;
            MPI_Type_create_subarray(2, gsizes, lsizes, start_ind, MPI_ORDER_C, MPI_DOUBLE, &filetype);
            MPI_Type_commit(&filetype);
          
            // Alterado para + 4
            int memsizes[2] = {myrows+4, mycols+4};
            start_ind[0] = 1;
            start_ind[1] = newid % 2;
            if (newid == 0 || newid == 1) {
                start_ind[0]--;
            }

            MPI_Datatype memtype;
            MPI_Type_create_subarray(2, memsizes, lsizes, start_ind, MPI_ORDER_C, MPI_DOUBLE, &memtype);
            MPI_Type_commit(&memtype);

            MPI_File fp;
            MPI_File_open(comm2D, "results_c.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
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
        // Alterado para +4
        MPI_Sendrecv(Vnew[myrows], mycols+4, MPI_DOUBLE, nbrtop, 0, Vnew[0] , mycols+4, MPI_DOUBLE, nbrbottom, 0, comm2D, MPI_STATUS_IGNORE);
        
        // comunicações sentido descendente
        // Alterado para usar 3a e antepenultima coluna (e +4)
        MPI_Sendrecv(Vnew[2], mycols+4, MPI_DOUBLE, nbrbottom, 1, Vnew[myrows+2] , mycols+4, MPI_DOUBLE, nbrtop, 1, comm2D, MPI_STATUS_IGNORE);
        
        // comunicações sentido para direita
        MPI_Sendrecv(&(Vnew[0][mycols]), 1, column, nbrright, 2, &(Vnew[0][0]), 1, column, nbrleft, 2, comm2D, MPI_STATUS_IGNORE);

        // comunicações sentido para esquerda
         // Alterado para usar 3a e antepenultima coluna
        MPI_Sendrecv(&(Vnew[0][2]), 1, column, nbrleft, 3, &(Vnew[0][mycols+2]), 1, column, nbrright, 3, comm2D, MPI_STATUS_IGNORE);
        

        // Outras comunicações
        // Comunicar 2a e penultima linhas
        MPI_Sendrecv(Vnew[myrows+1], mycols+4, MPI_DOUBLE, nbrtop, 6, Vnew[1], mycols+4, MPI_DOUBLE, nbrbottom, 6, comm2D, MPI_STATUS_IGNORE);

        // Comunicar 4a linha 
        MPI_Sendrecv(Vnew[3], mycols+4, MPI_DOUBLE, nbrbottom, 5, Vnew[myrows+3], mycols+4, MPI_DOUBLE, nbrtop, 5, comm2D, MPI_STATUS_IGNORE);
        
        // Comunicar 4a coluna 
        MPI_Sendrecv(&(Vnew[0][3]), 1, column, nbrleft, 7, &(Vnew[0][mycols+3]), 1, column, nbrright, 7, comm2D, MPI_STATUS_IGNORE);
        
        // Comunicar 2a e penultima colunas
        MPI_Sendrecv(&(Vnew[0][mycols+1]), 1, column, nbrright, 8, &(Vnew[0][1]), 1, column, nbrleft, 8, comm2D, MPI_STATUS_IGNORE);

        // Alterado para +4
        for (int i = 0; i < myrows + 4; i++)
        {
            for (int j = 0; j < mycols + 4; j++)
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