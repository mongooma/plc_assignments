#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


// Compile Code: mpicc -g -Wall main.c 
// Example Run Code: mpirun -np 4 ./a.out input.txt


int main(int argc, char **argv){

	int my_mpi_size;
	int my_mpi_rank; 

	MPI_Init( &argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	int * rank_arr = calloc(1, sizeof(int));
	int * tmp = calloc(1, sizeof(int));
	* tmp = -1;
	MPI_Request send_request;
	MPI_Request recv_request;
	MPI_Status	status;

	int arr[4] = {0, 1, 2, 3};
	int allocation = 1;


	MPI_Scatter(arr, allocation, MPI_INT, rank_arr, allocation, MPI_INT, 0, MPI_COMM_WORLD);

	//while(1){

		if(my_mpi_rank == 0){ //sender
			/*  int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
	            int tag, MPI_Comm comm, MPI_Request *request) */
			MPI_Isend(rank_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &send_request);
			//MPI_Wait(&send_request, &status);

		}else{

			while(1){
				MPI_Irecv(tmp, 1, MPI_INT, my_mpi_rank - 1, 0, MPI_COMM_WORLD, &recv_request);
				MPI_Wait(&recv_request, &status);
				if( * tmp != -1) break;
			}

			printf("rank %d recv %d \n", my_mpi_rank, * tmp);
			

			* rank_arr += * tmp;
			if(my_mpi_rank != (my_mpi_size -1)){
				MPI_Isend(rank_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &send_request);
				//MPI_Wait(&send_request, &status);
				printf("rank %d send %d \n", my_mpi_rank, * rank_arr);
			}
			

		}
		//if((my_mpi_rank == (my_mpi_size -1)) && * rank_arr != 0) break; 

		//sleep(1);
	//}

	MPI_Barrier(MPI_COMM_WORLD); 

	if(my_mpi_rank == (my_mpi_size -1)){
		printf("sum 0-3: %d\n", *rank_arr);
	}


	free(rank_arr);
	free(tmp);

	MPI_Finalize();

	return EXIT_SUCCESS;

	/* MPI_reduce could only perform sum on a single block and in parallel for other blocks*/

}