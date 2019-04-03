#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include"clcg4.h"

#include<mpi.h>
#include<pthread.h>

int main(int argc, char *argv[])
{

	int mpi_myrank;
    int mpi_commsize;
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

    // int ** recv = calloc(1, sizeof(int *));
    // int data[3][3] = {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}}; // on stack allocation
    int * data = calloc(3 * 3, sizeof(int));

    for(int i = 0; i < 9; i ++){
        data[i] = i;
    }

	int * recv = calloc(9, sizeof(int)); // row majored

    MPI_Request request_send;
    MPI_Request request_recv;

    // int flag;
    // MPI_Status status;




    /* int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request *request) */
    /* int MPI_Irecv(      void *buf, int count, MPI_Datatype datatype, int source, 
        int tag, MPI_Comm comm, MPI_Request *request) */

    if(mpi_myrank == 1){ // orders not matter

        // MPI_Send(*(data + 1), 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Isend(data + 3, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, &request_send);
        MPI_Wait(&request_send, MPI_STATUS_IGNORE); // block forever
    }
    
    if(mpi_myrank == 0){ // orders not matter

        MPI_Irecv(recv + 3, 3, MPI_INT, 1, 0, MPI_COMM_WORLD, &request_recv);
        MPI_Wait(&request_recv, MPI_STATUS_IGNORE); // block until some thing is received

    }




    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi_myrank == 0){
        printf("rank %d received:  %d, %d, %d\n", mpi_myrank, *(recv + 3), *(recv + 4), *(recv + 5));
    }

    MPI_Finalize();

    return 0;

}


