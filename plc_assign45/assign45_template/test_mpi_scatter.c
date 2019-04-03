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
    int * recv = calloc(3 * 3, sizeof(int)); // on stack allocation

    MPI_Barrier(MPI_COMM_WORLD);

	int data[4][3] = {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {5, 6, 7}}; // row majored

    // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //         void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    //         MPI_Comm comm)

    // MPI_Scatter(data, 3, MPI_INT, *(recv + 2), 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(data, 3, MPI_INT, recv + 3, 3, MPI_INT, 0, MPI_COMM_WORLD);
    // Scatter should be called from every rank ...

    MPI_Barrier(MPI_COMM_WORLD);

    printf("rank %d received:  %d, %d, %d\n", mpi_myrank, *(recv + 3), *(recv + 4), *(recv + 5));

    MPI_Finalize();

    return 0;

}