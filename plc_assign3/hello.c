// export PATH=/usr/local/mpich-3.2/bin:$PATH
// compile with: mpicc mpi-hello.c –o mpi-hello
// run: mpirun –np #pes mpi-hello
#include<stdio.h>
#include <mpi.h>

main(int argc, char *argv[])
{
	int npes, myrank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("From process %d out of %d, Hello World!\n",
		myrank, npes);
	MPI_Finalize();
}
