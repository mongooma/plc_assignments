#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>
#include "functions.h"


int main(int argc, char ** argv){
	/*
	if( argc != 2 ){
		printf("Not sufficient arguments, only %d found. \n", argc);
		exit(-1);
	}
	*/


	int my_mpi_size = -1; // total no. of ranks
	int my_mpi_rank = -1; // the current rank of the process

	MPI_Init( &argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	int bin1[4] = {0, 1, 2, 3};
	int bin2[4] = {0, 1, 2, 3};

	int * bin_rank_1 = calloc(4 / my_mpi_size, sizeof(int));
	int * bin_rank_2 = calloc(4 / my_mpi_size, sizeof(int));


	/*1. MPI rank 0 read in bit number, */

	//if(my_mpi_rank == 0){ // something wrong setting this line

	MPI_Scatter(bin1, 4 / my_mpi_size, MPI_INT, bin_rank_1, 4 / my_mpi_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(bin2, 4 / my_mpi_size, MPI_INT, bin_rank_2, 4 / my_mpi_size, MPI_INT, 0, MPI_COMM_WORLD);

	//} 


	for(int i = 0; i < 4 / my_mpi_size; i ++){
		printf("rank %d, bin1, %d, bin2, %d, \n", my_mpi_rank, 
				bin_rank_1[i], bin_rank_2[i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	free(bin_rank_1);
	free(bin_rank_2);

}
