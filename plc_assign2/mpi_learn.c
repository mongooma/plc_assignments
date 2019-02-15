#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>
#include "functions.h"


/* public cloud:
	

*/


int add(int my_mpi_rank, int * recv_arr){
	/* use some send and recv here*/

	/* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
               int source, int tag, MPI_Comm comm, MPI_Request *request) */
	

	int buf = 0 ;
	//MPI_Request request;
	MPI_Request request;

	// the following will not work, need to be explicitly written in order
	// if(my_mpi_rank == 0){
	// 	MPI_Isend(recv_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &status);
	// }

	// if(my_mpi_rank != 0){
	// 	MPI_Irecv(&buf, 1, MPI_INT, my_mpi_rank - 1, 0, MPI_COMM_WORLD, &status);
	// 	/* tags are pretty arbitrary in MPI msg */
	// 	/* MPI_ANY_TAG triggers error */


	// 	recv_arr[0] = buf + recv_arr[0];

	// 	if(my_mpi_rank != 3){
	// 		/*  int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
	// 	            int tag, MPI_Comm comm, MPI_Request *request)*/
	// 		MPI_Isend(recv_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &status);
	// 	}

	// }

	// method 1: use MPI_Send and MPI_Recv for blocking operation
	// for(int i = 0; i < 4; i ++){

	// 	if((my_mpi_rank == i) && (i == 0)){MPI_Send(recv_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD);}
		
	// 	if((my_mpi_rank == i) && (i != 0)){
	// 		MPI_Recv(&buf, 1, MPI_INT, my_mpi_rank - 1, 0, MPI_COMM_WORLD, &status); 
	// 		recv_arr[0] = buf + recv_arr[0];
			
	// 		if(i !=3 ){
	// 			MPI_Send(recv_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD);
	// 		}
	// 	}

	
	// }

	// method 2: use MPI_Isend and MPI_Irecv for non-blocking operation (return immediately)


	for(int i = 0; i < 4; i ++){

		if((my_mpi_rank == i) && (i == 0)){
			MPI_Isend(recv_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &request);
		}

		/* use barrier to have every ranks on the same level */
		MPI_Barrier(MPI_COMM_WORLD); /* don't put barrier in any if scope to cause blocking*/
		
		if((my_mpi_rank == i) && (i != 0)){
			MPI_Irecv(&buf, 1, MPI_INT, my_mpi_rank - 1, 0, MPI_COMM_WORLD, &request); 
			recv_arr[0] = buf + recv_arr[0];
		}

		MPI_Barrier(MPI_COMM_WORLD);
			
		if((my_mpi_rank == i) && (i < 3)){
			MPI_Isend(recv_arr, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &request);
		}

		MPI_Barrier(MPI_COMM_WORLD);

	}


	return 0;


}


/* mpicc -g -Wall mpi_learn.c */
/* mpiexec -np 4 ./a.out *//* using 4 ranks*/

int main(int argc, char ** argv){

	int my_mpi_size = -1;
	int my_mpi_rank = -1; /* somehow it's not hardcoded? */

	MPI_Init( &argc, &argv);
	/*  Open MPI accepts the C/C++ argc and
       argv arguments to main, but neither modifies, interprets, nor distributes them*/
  	
  	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
  	/*  Returns the size of the group associated with a communicator. */

  	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
  	/* Determines the rank of the calling process in the communicator.*/

  	
  	/*
  	int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm)
	*/

  	int arr[4] = {1, 2, 3, 4}; /* going to be distributed in the same order 1, 2, 3, 4*/
  	int * recv_arr = calloc(1, sizeof(int));
  	//if(my_mpi_rank == 0){ //something wrong if set this line? 
  	MPI_Scatter(arr, 1, MPI_INT, recv_arr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  	//}

  	printf("We are in rank %d, get %d from rank 0. \n", my_mpi_rank, recv_arr[0]);

  	MPI_Barrier(MPI_COMM_WORLD); /*  Synchronization between MPI processes in a group */

  	add(my_mpi_rank, recv_arr);

  	if(my_mpi_rank == 3){
	  	printf("The sum is: %d. \n", recv_arr[0]);
	  }

	/*  MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm) */
	int * arrs = calloc(4, sizeof(int));

	MPI_Gather(recv_arr, 1, MPI_INT, arrs, 1, MPI_INT, 0, MPI_COMM_WORLD); 

	MPI_Barrier(MPI_COMM_WORLD); /*  Synchronization between MPI processes in a group */
	/* Use barrier to get the right order of output */

	if(my_mpi_rank == 0){
		printf("We have gathered the modified recv_arrs: {%d...%d}. \n", arrs[0], arrs[3]);
	}

  	MPI_Finalize();

  	return 0;




}

int main_1(int argc, char ** argv){

	int my_mpi_size = -1;
	int my_mpi_rank = -1; /* somehow it's not hardcoded? */
	int * arr = calloc(4, sizeof(int));
	arr[0] = -1;


	MPI_Init( &argc, &argv);
	/*  Open MPI accepts the C/C++ argc and
       argv arguments to main, but neither modifies, interprets, nor distributes them*/
  	
  	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
  	/*  Returns the size of the group associated with a communicator. */

  	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
  	/* Determines the rank of the calling process in the communicator.*/

  	
  	/*
  	int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm)
	*/

  	
  	arr[0] = my_mpi_rank; // not using MPI_Scatter; each process has access with public data

  	printf("We are in rank %d, arr[0] is %d. \n", my_mpi_rank, arr[0]);


	MPI_Barrier(MPI_COMM_WORLD); /*  Synchronization between MPI processes in a group */
	/* Use barrier to get the right order of output */
  	MPI_Finalize();

  	free(arr);




}