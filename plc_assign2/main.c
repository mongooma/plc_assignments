#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>
#include "functions.h"

#define HEX_INPUT_SIZE 262144

/* variables for performance report:
	ranks: 2, 4, 8, 16, 32
	use_barrier: 0, 1

*/

// Compile Code: mpicc -g -Wall main.c 
// Example Run Code: mpirun -np 4 ./a.out input.txt


// Both input and output files will 524490 bytes in size. The two additional
//    characters are due to a newline characters between each input and at the
//    end of the file.



/*
 use 32 bit blocks and extend it to work in parallel for a 1M (e.g., 1,048,576) 
 bit CLA adder with 32 bit blocks using up to 32 MPI ranks on the class server, 
 mastiff.cs.rpi.edu. 
*/

/* 0. represent the 1,048,576 bit input numbers as 262,144 hex digits. and 
	and write to some file */

int main(int argc, char ** argv){
	/*
	if( argc != 2 ){
		printf("Not sufficient arguments, only %d found. \n", argc);
		exit(-1);
	}
	*/
	// Add 1 to array size because strings must be null terminated
	double start_time, end_time;

	int use_barrier = atoi(argv[argc-1]);
	//printf("Use use_barrier: %d", use_barrier);

	char * hex_input_a = calloc(HEX_INPUT_SIZE+1, sizeof(char)); //add 1 for '\0'
	char * hex_input_b = calloc(HEX_INPUT_SIZE+1, sizeof(char));

	//Integer array of inputs in binary form 4;
	int * bin1 = calloc((HEX_INPUT_SIZE) * 4, sizeof(int)); /* keep the fashion of assign 1, 0000 for the */
	int * bin2 = calloc((HEX_INPUT_SIZE) * 4, sizeof(int));


	int my_mpi_size = -1; // total no. of ranks
	int my_mpi_rank = -1; // the current rank of the process

	MPI_Init( &argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	int allocation = (HEX_INPUT_SIZE * 4) / my_mpi_size;
	int * bin_rank_1 = calloc(allocation, sizeof(int));
	int * bin_rank_2 = calloc(allocation, sizeof(int));
	int * sumi = calloc(allocation, sizeof(int));
	int * sumi_all = calloc(HEX_INPUT_SIZE * 4 , sizeof(int));

	/*1. MPI rank 0 read in bit number, */

	if(my_mpi_rank == 0){

		readInData(argv[1], hex_input_a, hex_input_b); /* not debugged */

		start_time = MPI_Wtime(); /* a time tick fashion*/

		/*2. MPI rank 0 convert hex to binary number, revert,  */

		convert_hex_2_bit(hex_input_a, hex_input_b, bin1, bin2, HEX_INPUT_SIZE); /* not debugged */

		revert_binary(bin1, bin2, HEX_INPUT_SIZE * 4);	

	}

	#ifdef DEBUG
	if(my_mpi_rank == 3){
		for(int i=0; i < 4; i ++){
			printf("bin1[%d]: %d \n", i, sumi[i]);
		}
	}
	#endif


	/*3. MPI Rank 0 distribute the input binary arrays to each rank in the correct order 
	where (for a 32 ranks conﬁguration) MPI rank 1 has bits 32,768 through 65,535 
	and MPI rank 2 has bits 65536 through 98304 and so on. */

	/* int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
        void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
        MPI_Comm comm)*/
	#ifdef DEBUG
		printf("rank %d, main: here00! \n", my_mpi_rank); /* three reached here */
	#endif 
	MPI_Scatter(bin1, allocation, MPI_INT, bin_rank_1, allocation, MPI_INT, 0, MPI_COMM_WORLD);
	#ifdef DEBUG
		printf("rank %d, main: here00! \n", my_mpi_rank); /* three reached here */
	#endif 
	MPI_Scatter(bin2, allocation, MPI_INT, bin_rank_2, allocation, MPI_INT, 0, MPI_COMM_WORLD);
	#ifdef DEBUG
		printf("rank %d, main: here01! \n", my_mpi_rank); /* three reached here */
	#endif 

	free(hex_input_a);
	free(hex_input_b);
	free(bin1);
	free(bin2);


	#ifdef DEBUG
		printf("rank %d, main: here1! \n", my_mpi_rank); /* three reached here */
	#endif 

	MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/



	/*execute algorithm -> see cla() */
	cla(use_barrier, my_mpi_rank, my_mpi_size, allocation, bin_rank_1, bin_rank_2, sumi); 

	#ifdef DEBUG
		printf("rank %d, main: here2! \n", my_mpi_rank); /* three reached here */
	#endif 


	/* for synchronization */
	MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/


	/* Have each rank send their part of the ﬁnal sumi solution to Rank 0. */

	/*  MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm) */	
	
	#ifdef DEBUG
	if(my_mpi_rank == 3){
		for(int i=0; i < 4; i ++){
			printf("sumi[%d]: %d \n", i, sumi[i]);
		}
	}
	#endif

	MPI_Gather(sumi, allocation, MPI_INT, sumi_all, allocation, MPI_INT, 0, MPI_COMM_WORLD); 
	

	/*  Rank 0 will then re-reverse the sum and output the ﬁnal result. */
	if( my_mpi_rank == 0){
		revert_hex_sum(sumi_all, HEX_INPUT_SIZE * 4);
		convert_bit_2_hex(sumi_all, HEX_INPUT_SIZE * 4);
		end_time = MPI_Wtime();
	}


	MPI_Finalize();


	if(my_mpi_rank == 0){ /*only monitor root process rank 0*/
		printf("time in seconds: %f, (use_barrier: %d, ranks: %s)\n", 
								(end_time - start_time), use_barrier, argv[0] );
	}

	free(bin_rank_1);
	free(bin_rank_2);
	free(sumi);
	free(sumi_all);


}