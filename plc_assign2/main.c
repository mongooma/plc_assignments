#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>
#include "functions.h"

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

/* 0. represent the 1,048,576 bit input numbers as 262,144 hex digits. and 
	and write to some file */

void readInData(){

	if( my_mpi_rank == 0 /* meaning the process in rank 0*/){

	#ifdef DEBUG
		printf("MPI Rank %d: Attempt to Read File Data \n", my_mpi_rank );
	#endif

		if( (my_input_file = fopen( argv[1], "r")) == NULL ){
			perror("Failed to open input data file: %s \n", argv[1]);
			exit(-1);
		}

		fscanf( my_input_file, "%s %s", hex_input_a, hex_input_b );
	
	#ifdef DEBUG
		printf("MPI Rank %d: Finished Reading and Write File Data \n", my_mpi_rank );
	#endif

		fclose( my_input_file );
	}

	MPI_Barrier(MPI_COMM_WORLD); /* collective operation*/

}


int main(int argc, char** argv){

	if( argc != 3 ){
		perror("Not sufficient arguments, only %d found \n", argc);
		exit(-1);
	}

	int my_mpi_size = ; // total no. of ranks
	int my_mpi_rank = ; // the current rank of the process

	MPI_Init( &argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD /* what is this? */, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	/*1. MPI rank 0 read in bit number, */

	readInData();

	/*2. NPI rank 0 convert bit to hex number */

	/*3. MPI Rank 0 distribute the input binary arrays to each rank in the correct order 
	where (for a 32 ranks conﬁguration) MPI rank 1 has bits 32,768 through 65,535 
	and MPI rank 2 has bits 65536 through 98304 and so on. */


	MPI_Scatter();

	/*execute algorithm -> see cla() */

	for(int i=0; i < my_mpi_size; i++){
		cla(); /* within it MPI_Barrier(MPI_COMM_WORLD) called from any rank will do?
				check man page */
	} 

	/* Have each rank send there part of the ﬁnal sumi solution to Rank 0. */

	MPI_Gather(); /* call this from any rank will do?*/

	MPI_Finalize();




}