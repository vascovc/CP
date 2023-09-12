// 97746 - Vasco Costa
// alinea D

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1e-6 // tolerancia
#define ITERMAX 500000 // criterio de paragem
#define NXMAX 500 
#define L 1.0 // Tamanho do lado do quadrado, do dominio

double f(double x, double y){
    //return (x + y);
    return 7*sin(2*M_PI*x)*cos(3*M_PI*x)*sin(2*M_PI*y)*cos(3*M_PI*y);
}

int main(int argc, char *argv[])
{
    int nprocs; // Numero de processos
    int myid; // Rank  de cada processo
    int nx, ny; // Tamanho da matriz para ambas as coordenadas
    double t0;   // Tempo inicial
    
    // Definir o MPI
    MPI_Init(&argc, &argv); // Inicializar o MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
        printf("Introduza o numero de pontos (max %d, 0 para sair)", NXMAX);
        scanf(" %d", &nx);
        //nx = 100;
        t0 = MPI_Wtime();
    }

    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    ny = nx;

    if(nx == 0){ // O utilizador quer sair
        MPI_Finalize();
        return 0;
    }
    else if (nx < 0 || nx > NXMAX) // invalido ou demasiado grande
    {
        MPI_Finalize();
        return 1;
    }

    // Definicao e inicializacao de variaveis para a criacao do comunicador cartesiano
    int ndims = 2; // Numero de dimensoes
    int dims[2] = {(int)(nprocs / 2), 2}; // Numero de linhas, numero de colunas
    
    // alterado
    int periodic[2] = {1, 1}; 
    // fim alteracao
    
    MPI_Comm comm2D;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, 1, &comm2D); // O 1 e para reordenar os ranks. O &comm2D e para guardar o novo comunicador
    if(comm2D == MPI_COMM_NULL){ // Se o comunicador for nulo, entao o processo nao pertence ao comunicador
        MPI_Finalize();
        return 0;
    }
    // Definir o novo id, isto e, o novo rank
    int newid;
    MPI_Comm_rank(comm2D, &newid); // O newid vai ser o rank do processo no novo comunicador

    // Definir os vizinhos
    int nbrbottom, nbrtop, nbrright, nbrleft;
    // Cada processo tem de saber o rank dos seus vizinhos, o de cima e o de baixo
    MPI_Cart_shift(comm2D, 0, 1, &nbrtop, &nbrbottom); // Fazemos o shift no comm2D para descobrir vizinho de cima e baixo
    MPI_Cart_shift(comm2D, 1, 1, &nbrleft, &nbrright); // Fazemos o shift no comm2D para descobrir vizinho da direita e da esquerda

    // verificar os vizinhos - debug
    printf("newid = %d   nbrtop = %d     nbrbottom = %d     nbrleft = %d    nbrright = %d\n", newid, nbrtop, nbrbottom, nbrleft, nbrright);
    MPI_Barrier(comm2D);

    int firstrow, nrows;
    int firstcol, ncols; // Para a versao 2D

    // efeciencia
    int dim0 = dims[0];
    // Atualizar o numero de processos
    nprocs = dim0 * dims[1];

    // O processo 0, porque este existe sempre, e quem cria estes arrays
    if (newid == 0)
    {   
        int listfirstrow[nprocs];
        int listnrows[nprocs];

        int listfirstcol[nprocs];
        int listncols[nprocs];

        int nrowsP = (int)((double)(ny-2)/(double)dim0 + 0.5);  

        // Ciclo for para percorrer todos os processos. O i e o rank do processo
        /*
        for (int i = 0; i < dim0; i++){
            // Para cada processo indicamos
            listfirstrow[2*i] = 1 + i *  nrowsP;
            listnrows[2*i] = nrowsP;
            listfirstrow[2*i+1] = 1 + i *  nrowsP;
            listnrows[2*i+1] = nrowsP;
        }
        // altera o numero de linhas do penultimo e do ultimo
        listnrows[nprocs-2] = ny - 2 - (dim0 - 1) * nrowsP;
        listnrows[nprocs-1] = ny - 2 - (dim0 - 1) * nrowsP;
        // Colunas
        int ncols_temp = (int)((nx-2)/2);
        for (int i = 0; i < dim0; i++){
            //pares
            listfirstcol[2*i] = 1;
            listncols[2*i] = ncols_temp;
            //impares
            listfirstcol[2*i+1] = ncols_temp + 1;
            listncols[2*i+1] = nx - 2 - ncols_temp;
        }
        */

        // alteracoes
        for (int i = 0; i < dim0; i++)
        {
            listfirstrow[2*i] = i * nrowsP;
            listnrows[2*i] = nrowsP+1; 
            listfirstrow[2*i+1] = i * nrowsP;
            listnrows[2*i+1] = nrowsP+1;
        }

        listfirstrow[nprocs - 1] = 1 + (dim0 - 1) * nrowsP; 
        listnrows[nprocs - 1] = ny - (dim0 - 1) * nrowsP;

        listfirstrow[nprocs - 2] = 1 + (dim0 - 1) * nrowsP;  
        listnrows[nprocs - 2] = ny - (dim0 - 1) * nrowsP;

        int ncols_temp = (int)((nx-2)/2);
        for (int i = 0; i < dim0; i++)
        {
            listfirstcol[2*i] = 0;
            listncols[2*i] = ncols_temp + 1;
            listfirstcol[2*i+1] = ncols_temp + 1;
            listncols[2*i+1] = nx - 1 - ncols_temp;
        }
        //fim alteracoes

        //linhas
        MPI_Scatter(listfirstrow, 1, MPI_INT, &firstrow, 1, MPI_INT, newid, comm2D);
        MPI_Scatter(listnrows, 1, MPI_INT, &nrows, 1, MPI_INT, newid, comm2D);
        //colunas
        MPI_Scatter(listfirstcol, 1, MPI_INT, &firstcol, 1, MPI_INT, newid, comm2D);
        MPI_Scatter(listncols, 1, MPI_INT, &ncols, 1, MPI_INT, newid, comm2D);

    }
    else{
        // Se nao for o processo 0, entao vai receber as linhas e colunas que lhe sao atribuidas
        // Aqui so se vai receber, apontador NULL (MPI_BOTTOM) que nao envia nada. e sempre o 0 que envia. Ele esta aqui para os outros processos saberem que nao sao eles que enviam. Eles tem o papel de receber
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstrow, 1, MPI_INT, 0, comm2D);
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &nrows, 1, MPI_INT, 0, comm2D);

        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstcol, 1, MPI_INT, 0, comm2D);
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &ncols, 1, MPI_INT, 0, comm2D);
    }

    //debug
    printf("newid = %d   firstrow = %d     lastrow = %d     firstcol = %d       lastcol = %d\n", newid, firstrow, firstrow + nrows - 1, firstcol, firstcol + ncols - 1);
    MPI_Barrier(comm2D);

    // Declarar matrizes com declaracaoo dinamica
    double (*Vnew)[ncols + 2]; // Numero de linhas da matriz Vnem. O +2 e por causa das fantasma
    double (*Vold)[ncols + 2];
    double (*myf)[ncols + 2];
    Vnew = calloc(nrows + 2, sizeof(*Vnew)); // O numero de linhas tem de ser igual ao numero de linhas +2 -> "fantastma"                                 
    Vold = calloc(nrows + 2, sizeof(*Vold));
    myf = calloc(nrows + 2, sizeof(*myf));

    double h = ((double)2 * L) / ((double) nx);
    for (int j = 1; j <= ncols ; j++){
        for (int i = 1; i <= nrows ; i++)
        {
            myf[i][j] = f(-L + (firstcol + j - 1) * h, -L + (firstrow + i - 1) * h);
        }
        
    }

    // condicoes fronteira originais
    /*
                // Preencher as condicoes fronteira. Esses valores nao sao alterados mas necessitam de ser inicializados
                // Para o processo de cima/linha de cima que agora esta no processo 0 e 1
                if(newid == 0 || newid == 1){
                    // Percorrer elementos da linha de cima
                    for(int i = 0; i < (ncols + 2); i++){ // Ou comecar em i = 1 e ter i < ncols
                        // Preeenche o valor da condicao fronteira, elemento a elemento
                        Vnew[0][i] = 0;
                        Vold[0][i] = 0;
                    }
                }
                // Para o processo de baixo/linha de baixo
                if(newid == (nprocs - 1) || newid == (nprocs - 2)){
                    // Percorrer elementos da linha de baixo
                    for(int i = 0; i < (ncols + 2); i++){
                        // Preeenche o valor da condicaoo fronteira, elemento a elemento
                        Vnew[nrows + 1][i] = 0;
                        Vold[nrows + 1][i] = 0;
                    }
                }

                // Se der 0 e par, se der 1 e impar. Fazemos isto porque as condicoes fronteira sao diferentes para os processo par e impar
                if(newid % 2 == 0){
                    for(int i = 0; i < (nrows + 2); i++){
                        // Preeenche o valor da condicao fronteira da esquerda, elemento a elemento
                        Vnew[i][0] = 0;
                        Vold[i][0] = 0;
                    }
                }
                else{
                    for(int i = 0; i < (nrows + 2); i++){
                        // Preeenche o valor da condicao fronteira da direita, elemento a elemento
                        Vnew[i][ncols + 1] = 0;
                        Vold[i][ncols + 1] = 0;
                    }

                }
    */
    MPI_Datatype column;
    MPI_Type_vector((nrows + 2), 1, (ncols + 2) , MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    
    // Iterar por cada dominio - Iterar o momento Jacobi
    for (int iter = 0; iter < ITERMAX; iter++){
        // Definir um array para somas
        double sums[2] = {0.0,0.0};
        double global_sums[2];

        /*
        for (int j = 1; j < ncols + 1 ; j++)
        {
            for (int i = 1; i < nrows + 1; i++)
            {
                Vnew[i][j] = (Vold[i+1][j] + Vold[i-1][j] + Vold[i][j+1] + Vold[i][j-1] + h * h  * myf[i][j]) / 4.0;
                sums[0] += (Vnew[i][j] - Vold[i][j]) * (Vnew[i][j] - Vold[i][j]);
                sums[1] += Vnew[i][j] * Vnew[i][j];
            }
            
        }
        */
        //alteracoes
        // calulo para os pares
        for (int i = 1; i < nrows + 1; i++){
            for (int j = 1; j < ncols + 1 ; j++){
                if (((firstcol + j - 1)  + (firstrow + i - 1)) %2 == 0){ // par
                    Vnew[i][j] = (Vnew[i+1][j] + Vnew[i-1][j] + Vnew[i][j+1] + Vnew[i][j-1]  + h * h  * myf[i][j]) / 4.0;
                    sums[0] += (Vnew[i][j] - Vold[i][j]) * (Vnew[i][j] - Vold[i][j]);
                    sums[1] += Vnew[i][j] * Vnew[i][j];
                }
            }  
        }
        // comunicar aos vizinhos
        MPI_Sendrecv(&Vnew[1][1], ncols, MPI_DOUBLE, nbrtop, 4, &Vnew[nrows+1][1], ncols, MPI_DOUBLE, nbrtop, 4, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&Vnew[nrows][1], ncols, MPI_DOUBLE, nbrbottom, 5, &Vnew[0][1], ncols, MPI_DOUBLE, nbrtop, 5, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&Vnew[1][1], 1, column, nbrleft, 6, &Vnew[1][ncols+1], 1, column, nbrright, 6, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&Vnew[1][ncols], 1, column, nbrright, 7, &Vnew[1][0], 1, column, nbrleft, 7, comm2D, MPI_STATUS_IGNORE);

        // impar
        for (int i = 1; i < nrows + 1; i++){
            for (int j = 1; j < ncols + 1 ; j++){
                if (((firstcol + j - 1)  + (firstrow + i - 1)) %2 == 1){ // impar
                    Vnew[i][j] = (Vnew[i-1][j] + Vnew[i][j-1] + Vnew[i][j+1] + Vnew[i+1][j] + h * h * myf[i][j]) / 4.0 ;
                    sums[0] += (Vnew[i][j]-Vold[i][j])*(Vnew[i][j]-Vold[i][j]);
                    sums[1] += Vnew[i][j]*Vnew[i][j];
                }
            }
        }

        // comunicar aos vizinhos
        MPI_Sendrecv(&Vnew[1][1], ncols, MPI_DOUBLE, nbrtop, 8, &Vnew[nrows+1][1], ncols, MPI_DOUBLE, nbrbottom, 8, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&Vnew[nrows][1], ncols, MPI_DOUBLE, nbrbottom, 9, &Vnew[0][1], ncols, MPI_DOUBLE, nbrtop, 9, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&Vnew[1][1], 1, column, nbrleft, 10, &Vnew[1][ncols+1], 1, column, nbrright, 10, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&Vnew[1][ncols], 1, column, nbrright, 11, &Vnew[1][0], 1, column, nbrleft, 11, comm2D, MPI_STATUS_IGNORE);
        // fim alteracoes

        MPI_Allreduce(sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, comm2D);
        
        if (sqrt(global_sums[0]/global_sums[1]) < TOL)
        {
            if(newid==0){
				printf("Calculo terminado com %d iteracoes\n",iter);
				printf("Tempo de execucao = %f\n",MPI_Wtime()-t0);
                t0 = MPI_Wtime();
			}
            int gsizes[2] = {ny, nx};
            int lsizes[2] = {nrows, ncols}; 
            int start_ind[2] = {firstrow, firstcol};
            // ja nao e necessario ter em conta as colunas fantasma   
            /*
            int lsizes[2] = {nrows, ncols + 1};
            int start_ind[2] = {firstrow, firstcol - 1 + (newid % 2)};
            if (newid == 0 || newid == 1) {
                lsizes[0]++;
                start_ind[0]--;
            }
            if (newid == nprocs-2 || newid == nprocs-1) {
                lsizes[0]++;
            }
            */
            MPI_Datatype filetype;
            MPI_Type_create_subarray(2, gsizes, lsizes, start_ind, MPI_ORDER_C, MPI_DOUBLE, &filetype);
            MPI_Type_commit(&filetype);
          
            int memsizes[2] = {nrows+2, ncols+2};
            start_ind[0] = 1;
            start_ind[1] = newid % 2;
            if (newid == 0 || newid == 1) {
                start_ind[0]--;
            }

            MPI_Datatype memtype;
            MPI_Type_create_subarray(2, memsizes, lsizes, start_ind, MPI_ORDER_C, MPI_DOUBLE, &memtype);
            MPI_Type_commit(&memtype);

            MPI_File fp;
            MPI_File_open(comm2D, "results_2D_d.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
            MPI_File_set_view(fp, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
            
            MPI_File_write_all(fp, Vnew, 1, memtype, MPI_STATUS_IGNORE);
            MPI_File_close(&fp);
            if(newid==0){
                printf("Tempo de escrita = %f\n",MPI_Wtime()-t0);
			}

            MPI_Type_free(&filetype);
            MPI_Type_free(&memtype);
            break;
        }

        // trocar linhas fantasmas
        MPI_Sendrecv(Vnew[nrows], ncols+2, MPI_DOUBLE, nbrbottom, 0, Vnew[0] , ncols+2, MPI_DOUBLE, nbrtop, 0, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(Vnew[1], ncols+2, MPI_DOUBLE, nbrtop, 1, Vnew[nrows+1] , ncols+2, MPI_DOUBLE, nbrbottom, 1, comm2D, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&(Vnew[0][ncols]), 1, column, nbrright, 2, &(Vnew[0][0]), 1, column, nbrleft, 2, comm2D, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&(Vnew[0][1]), 1, column, nbrleft, 3, &(Vnew[0][ncols+1]), 1, column, nbrright, 3, comm2D, MPI_STATUS_IGNORE);
        
        // copiar o Vnew para o Vold
        for (int i = 0; i < nrows + 2; i++)
        {
            for (int j = 0; j < ncols + 2; j++)
            {
                Vold[i][j] = Vnew[i][j];
            }
            
        }
        
    }

    MPI_Type_free(&column);
    free(Vold);
    free(Vnew);
    free(myf);
    MPI_Comm_free(&comm2D);
    MPI_Finalize();

    return 0;
}