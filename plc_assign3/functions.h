#include <mpi.h>


int MPI_P2P_reduce(const unsigned long long *sendbuf, unsigned long long *recvbuf, int count,
                       MPI_Datatype datatype, MPI_Op op, int root,
                       MPI_Comm comm);