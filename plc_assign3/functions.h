#include <mpi.h>


int MPI_P2P_reduce(const int *sendbuf, int *recvbuf, int count,
                       MPI_Datatype datatype, MPI_Op op, int root,
                       MPI_Comm comm);