#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TOOL 1e-5
#define ITER_MAX 1000
#define NXMAX 500
#define L 1.0

double f(double x, double y){
	return x+y;
}

int main(int argc, char *argv[]){

	int nprocs;
	int myid;
	int nx, ny;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if(myid == 0){
		printf("Numero de pontos (max %d, 0 para sair): ",NXMAX);
		scanf(" %d", &nx);
	}

	MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ny=nx;

	if(nx == 0) {
		MPI_Finalize();
		return 0;
	}
	else if (nx > NXMAX){
		printf("Numero de pontos superior ao permitido\n");
		MPI_Finalize();
		return 1;
	}

	int ndims = 2;
	int dims[2] = {(int)(nprocs/2),2};
	int periodic[2] = {0,0};
	MPI_Comm comm2d;
	int newid,nbrbottom,nbrtop,nbrleft,nbrright;

	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic,1, &comm2d);
	MPI_Comm_rank(comm2d, &newid);

	MPI_Cart_shift(comm2d, 0, 1, &nbrtop, &nbrbottom);
    MPI_Cart_shift(comm2d, 1, 1, &nbrleft, &nbrright);

	printf("newid=%d    nbrtop=%d     nbrbottom=%d     nbrleft=%d     nbrright=%d\n",newid,nbrtop,nbrbottom,nbrleft,nbrright);
    MPI_Barrier(comm2d);
	int firstrow,nrows;
    int firstcol,ncols;

    if(comm2d == MPI_COMM_NULL){
        MPI_Finalize();
        return 0;
    }
    int dims0=dims[0];
    nprocs = dims[0]*dims[1];

	if(newid == 0){
		int listfirstrow[nprocs];
		int listnrows[nprocs];
        int listfirstcol[nprocs];
        int listncols[nprocs];
		
		int nrows_p = (int)(((double)(ny-2))/((double)dims0)+0.5);

		for(int i=0;i<dims0-1;i++){
			listfirstrow[2*i] = 1+i*nrows_p;
			listnrows[2*i] = nrows_p;
            listfirstrow[2*i+1] = 1+i*nrows_p;
			listnrows[2*i+1] = nrows_p;
		}
		listfirstrow[nprocs-1] = 1+(dims0-1)*nrows_p;
		listnrows[nprocs-1] = ny-2-(dims0-1)*nrows_p;
        listfirstrow[nprocs-2] = 1+(dims0-1)*nrows_p;
		listnrows[nprocs-2] = ny-2-(dims0-1)*nrows_p;

        for(int i=0;i<dims0;i++){
            listfirstcol[2*i] = 1;
            listncols[2*i] = (int)(nx/2)-1;
            listfirstcol[2*i+1] = (int)(nx/2);
            listncols[2*i+1] = ny-1-(int)(nx/2);
        }

		MPI_Scatter(listfirstrow, 1, MPI_INT, &firstrow, 1, MPI_INT, 0, comm2d);
		MPI_Scatter(listnrows, 1, MPI_INT, &nrows, 1, MPI_INT, 0, comm2d);
        MPI_Scatter(listfirstcol, 1, MPI_INT, &firstcol, 1, MPI_INT, 0, comm2d);
		MPI_Scatter(listncols, 1, MPI_INT, &ncols, 1, MPI_INT, 0, comm2d);
	}else{
		//envia o 0 e aqui fica o 0 que e para saberem que n e ele que esta a enviar
		MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstrow, 1, MPI_INT, 0, comm2d);
		MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &nrows, 1, MPI_INT, 0, comm2d);
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstcol, 1, MPI_INT, 0, comm2d);
		MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &ncols, 1, MPI_INT, 0, comm2d);
}

	printf("newid=%d  firstrow=%d  lastrow=%d  firstcol=%d lastcol=%d\n",newid,firstrow,firstrow+nrows-1,firstcol,firstcol+ncols-1);
    MPI_Barrier(comm2d);

	double(*Vnew)[ncols+2], (*Vold)[ncols+2], (*myf)[ncols+2];
	Vnew = calloc(nrows+2,sizeof(*Vnew));
	Vold = calloc(nrows+2,sizeof(*Vold));
	myf = calloc(nrows+2,sizeof(*myf));

	double h = L/((double)(nx-1));
	//inicializar matriz da funcao f
	for(int i=1;i<=nrows;i++){
		for(int j=1;j<=ncols;j++){
			myf[i][j] = f((double)(-L/2+(firstcol+j-1)*h),(double)(-L/2-(firstrow+i-1)*h));
		}
	}

	// inicializar condicoes fronteira
	if (newid == 0 || newid == 1){
		for(int i=0;i<ncols+2;i++){
			Vnew[0][i] = 0;
			Vold[0][i] = 0;
		}
	}
	if (newid == nprocs-1 || newid == nprocs-2){
		for(int i=0;i<nx;i++){
			Vnew[nrows+1][i] = 0;
			Vold[nrows+1][i] = 0;
		}
	}
    if(newid % 2==0){
        for (int i = 0; i < nrows+2; i++)
        {
            Vnew[i][0] = 0;
            Vold[i][0] = 0;
        }
    } else {
        for (int i = 0; i < nrows+2; i++)
        {
            Vnew[i][ncols+1] = 0;
            Vold[i][ncols+1] = 0;
        }
    }

    MPI_Datatype column;
    MPI_Type_vector((nrows+2), 1, (ncols+2), MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

	for(int iter=0;iter<ITER_MAX;iter++){
		double sums[2];
		for(int i=1;i<=nrows;i++){
			for(int j=1;j<=ncols-1;j++){
				Vnew[i][j] = 0.25*(Vold[i-1][j]+Vold[i+1][j]+Vold[i][j-1]+Vold[i][j+1]-h*h*myf[i][j]);
				sums[0] += (Vnew[i][j]-Vold[i][j])*(Vnew[i][j]-Vold[i][j]);
				sums[1] += Vnew[i][j]*Vnew[i][j];
			}
		}
		double global_sums[2];
		// troca de informacao entre processos
		
		
		MPI_Allreduce(sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, comm2d);
		if(global_sums[0]/global_sums[1] < TOOL){
			//escrever
			/*
			printf("p%d\n", newid);
			for(int i=0;i<nrows+2;i++){
				for(int j=0;j<ncols+2;j++){
					printf("%f ", Vnew[i][j]);
				}
				printf("\n");
			}
			*/
			break;
		}
		
		// troca de informacao entre processos
		/*
		MPI_Allreduce(sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, comm2d);
		if(global_sums[0]/global_sums[1] < TOOL){
			//escrever
			if(newid==0){
				printf("Calculo terminado com %d iteracoes\n",iter);
				//printf("Tempo de execucao = %f\n",MPI_Wtime()-t0);
			}
			int gsizes[2] = {nx,ny};
			int lsizes[2] = {nrows,ncols+1};
			int start_ind[2] = {firstrow,firstcol-1+newid%2};
			if(newid==0 || newid == 1){
				lsizes[0]++;
				start_ind[0]--;
			}
			if(newid==nprocs-1 || newid == nprocs-2){
				lsizes[0]++;
			}
			MPI_Datatype filetype;
			MPI_Type_create_subarray(2,gsizes,lsizes,start_ind,MPI_ORDER_C,MPI_DOUBLE,&filetype);
			MPI_Type_commit(&filetype);

			int memsizes[2] = {nrows+2,ncols+2};
			start_ind[0] = 1;
			start_ind[1] = newid%2;
			if(newid==0 || newid == 1){
				start_ind[0]=0;
			}
			MPI_Datatype memtype;
			MPI_Type_create_subarray(2,memsizes,lsizes,start_ind,MPI_ORDER_C,MPI_DOUBLE,&memtype);
			MPI_Type_commit(&memtype);

			MPI_File fp;
			MPI_File_open(comm2d,"results_2D.bin",MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fp);
			MPI_File_set_view(fp,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);

			MPI_File_write_all(fp,Vnew,1,memtype,MPI_STATUS_IGNORE);

			MPI_File_close(&fp);

			MPI_Type_free(&filetype);
			MPI_Type_free(&memtype);


			break;
		}
*/

		//sentido ascendente
		MPI_Sendrecv(Vnew[1], (ncols + 2), MPI_DOUBLE, nbrtop, 0, Vnew[nrows + 1], (ncols + 2), MPI_DOUBLE, nbrbottom, 0, comm2d, MPI_STATUS_IGNORE); // Sentido ascendente
        MPI_Sendrecv(Vnew[nrows], (ncols + 2), MPI_DOUBLE, nbrbottom, 1, Vnew[0], (ncols + 2), MPI_DOUBLE, nbrtop, 1, comm2d, MPI_STATUS_IGNORE); // Sentido descendente
        
        MPI_Sendrecv(&(Vnew[0][ncols]), 1, column, nbrright, 2, &(Vnew[0][0]), 1, column, nbrleft, 2, comm2d, MPI_STATUS_IGNORE); // Sentido para a direita
        MPI_Sendrecv(&(Vnew[0][1]), 1, column, nbrleft, 3, &(Vnew[0][ncols + 1]), 1, column, nbrright, 3, comm2d, MPI_STATUS_IGNORE); // Sentido para a esquerda

		//atualizar Vold
		for(int i=0;i<nrows+2;i++){
			for(int j=0;j<ncols+2;j++){
				Vold[i][j] = Vnew[i][j];
			}
		}


	}

    MPI_Type_free(&column);
	free(Vnew), free(Vold), free(myf);
	MPI_Finalize();
	return 0;
}
