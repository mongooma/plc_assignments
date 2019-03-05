#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


// Compile Code: mpicc -g -Wall main.c 
// Example Run Code: mpirun -np 4 ./a.out input.txt


int main(int argc, char **argv){

	int my_mpi_size;
	int my_mpi_rank; 

	MPI_Init( &argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	int arr[8] = {0, 1, 2, 3};
	int allocation = 1;
	int * rank_arr = calloc(1, sizeof(int));
	int sum;
	MPI_Scatter(arr, allocation, MPI_INT, rank_arr, allocation, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Reduce(rank_arr, &sum, allocation,
                       MPI_INT, MPI_SUM, 0,
                       MPI_COMM_WORLD);

	if(my_mpi_rank == 0){
		printf("sum 0-8: %d\n", sum);
	}


	free(rank_arr);

	MPI_Finalize();

	return EXIT_SUCCESS;

	/* MPI_reduce could only perform sum on a single block and in parallel for other blocks*/

}

