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
    int recv[3][3] = {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}}; // on stack allocation

    printf("recv[2][0]: %d\n", *(recv + 2)[0]);

    MPI_Barrier(MPI_COMM_WORLD);

	int data[4][3] = {{0}}; // row majored
    MPI_Scatter(data, 3, MPI_INT, *(recv + 2), 3, MPI_INT, 0, MPI_COMM_WORLD);
    // Scatter should be called from every rank ...

    MPI_Barrier(MPI_COMM_WORLD);

    printf("rank %d received:  %d, %d, %d\n", mpi_myrank, recv[2][0], recv[2][1], recv[2][2]);

    MPI_Finalize();

    return 0;

}