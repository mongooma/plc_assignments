#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>
#include "functions.h"

#define HEX_INPUT_SIZE 262144
#define use_barrier 1

/* variables for performance report:
	ranks: 2, 4, 8, 16, 32
	use_barrier: 0, 1

*/

// Compile Code: mpicc -g -Wall mpi-cla-io.c -o mpi-cla-io
// Example Run Code: mpirun -np 4 ./mpi-cla-io test_input_1.txt test_output_1.txt


// Both input and output files will 524490 bytes in size. The two additional
//    characters are due to a newline characters between each input and at the
//    end of the file.



/*
 use 32 bit blocks and extend it to work in parallel for a 1M (e.g., 1,048,576) 
 bit CLA adder with 32 bit blocks using up to 32 MPI ranks on the class server, 
 mastiff.cs.rpi.edu. 
*/

FILE *my_input_file=NULL;

// Add 1 to array size because strings must be null terminated
char hex_input_a[HEX_INPUT_SIZE+1]={0};
char hex_input_b[HEX_INPUT_SIZE+1]={0};

//Integer array of inputs in binary form 
bin1 = calloc((HEX_INPUT_SIZE+1) * 4, sizeof(int)); /* keep the fashion of assign 1, 0000 for the */
bin2 = calloc((HEX_INPUT_SIZE+1) * 4, sizeof(int));


/* 0. represent the 1,048,576 bit input numbers as 262,144 hex digits. and 
	and write to some file */

int main(int argc, char** argv){

	if( argc != 3 ){
		perror("Not sufficient arguments, only %d found \n", argc);
		exit(-1);
	}

	int my_mpi_size = -1; // total no. of ranks
	int my_mpi_rank = -1; // the current rank of the process

	double start_time = MPI_Wtime(); /* a time tick fashion*/

	MPI_Init( &argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	int allocation = ((HEX_INPUT_SIZE + 1) * 4) / my_mpi_size;
	int * bin_rank_1 = calloc(allocation, sizeof(int));
	int * bin_rank_2 = calloc(allocation, sizeof(int));
	int * sumi = calloc(allocation, sizeof(int));
	int * sumi_all = calloc((HEX_INPUT_SIZE+1) * 4 , sizeof(int));

	/*1. MPI rank 0 read in bit number, */

	if(my_mpi_rank == 0){

		readInData(argv[1], hex_input_a, hex_input_b); /* not debugged */

		/*2. MPI rank 0 convert hex to binary number, revert,  */

		convert_hex_2_bit(hex_input_a, hex_input_b, bin1, bin2, HEX_INPUT_SIZE); /* not debugged */

		revert_binary(bin1, bin2, (HEX_INPUT_SIZE + 1) * 4);	

		/*3. MPI Rank 0 distribute the input binary arrays to each rank in the correct order 
		where (for a 32 ranks conﬁguration) MPI rank 1 has bits 32,768 through 65,535 
		and MPI rank 2 has bits 65536 through 98304 and so on. */

		/* int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm)*/
		MPI_Scatter(bin1, allocation, MPI_INT, bin_rank_1, allocation, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(bin2, allocation, MPI_INT, bin_rank_2, allocation, MPI_INT, 0, MPI_COMM_WORLD);

	}


	/*execute algorithm -> see cla() */
	cla(use_barrier, my_mpi_rank, my_mpi_size, allocation, bin_rank_1, bin_rank_2); 

	/* for synchronization */
	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}


	/* Have each rank send their part of the ﬁnal sumi solution to Rank 0. */

	/*  MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm) */	
	MPI_Gather(sumi, allocation, MPI_INT, sumi_all, allocation, MPI_INT, 0, MPI_COMM_WORLD); 

	/*  Rank 0 will then re-reverse the sum and output the ﬁnal result. */
	if( my_mpi_rank == 0){
		revert_hex_sum(sumi_all, (HEX_INPUT_SIZE + 1) * 4);
		convert_bit_2_hex(sumi_all, (HEX_INPUT_SIZE + 1) * 4);

	}

	double end_time = MPI_Wtime();

	MPI_Finalize();


	printf("time in seconds: %d, (use_barrier: %d, ranks: %d)", 
								end_time, use_barrier, argv[0]);


	free(hex_input_a);
	free(hex_input_b);
	free(bin1);
	free(bin2);
	free(bin_rank_1);
	free(bin_rank_2);
	free(sumi);
	free(sumi_all);


}