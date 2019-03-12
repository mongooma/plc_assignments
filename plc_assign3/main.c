#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "functions.h"

// #define BGQ 1 // when running BG/Q, comment out when running on mastiff
#ifdef BGQ 
#include <hwi/include/bqc/A2_inlines.h> 
#else 
#define GetTimeBase MPI_Wtime 
#endif


/*
	For this assignment, your are to develop a point-2-point message version of 
	the MPI_reduce operation and compare the performance of that version you develop 
	to the collective version, MPI_reduce across a variable number of MPI rank 
	configurations.

*/

/*
	Memory tips:
	1 node run will than have 64 arrays of 16,777,216 entries each. 
	However, since each array entry is an 8 byte quantity, 
	it will consume 134,217,728 bytes per rank which is about 1/2 the 
	available memory on each rank.

*/

/*
	Experiment:

	You will conduct a strong scaling experiment using “AMOS”, the CCI Blue Gene/Q system. Here, you will have 64 MPI rank per Blue Gene/Q node.
	• Run 1 billion entry 64 MPI ranks ( 1 BG/Q nodes )
	• Run 1 billion entry 128 MPI ranks ( 2 BG/Q nodes )
	• Run 1 billion entry 256 MPI ranks ( 4 BG/Q nodes )
	• Run 1 billion entry 512 MPI ranks ( 8 BG/Q nodes )
	• Run 1 billion entry 1024 MPI ranks ( 16 BG/Q nodes )
	• Run 1 billion entry 2048 MPI ranks ( 32 BG/Q nodes )
	• Run 1 billion entry 4096 MPI ranks ( 64 BG/Q nodes )
	• Run 1 billion entry 8192 MPI ranks ( 128 BG/Q nodes )

*/

/* For different ranks, may want to run with an environment number setting:
	env RANK_no=... 

	Since I haven't found solutions to "extracting the no. of ranks within the communicator during running"
	*/


int main(int argc, char ** argv){

	/* Use the Blue Gene/Q’s “GetTimeBase()” function to measure the number of 
		cycles this point-2-point reduce operation took as well as the MPI 
		collection MPI reduce
	*/
	double time_in_secs = 0; 
	double time_in_secs_mpi = 0; 
	double processor_frequency = 1600000000.0; 
	unsigned long long start_cycles=0; 
	unsigned long long end_cycles=0; 
	int my_mpi_size;
	int my_mpi_rank; 
	MPI_Init( &argc, &argv); /* also take params from regular command line input */
	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	unsigned long long sum = 0;
	unsigned long long sum_mpi = 0;

	/* demo */
	//unsigned long long arr[8] = {0, 1, 2, 3, 4, 5, 6, 7};
	//int allocation = 2;
	/***/

	// MPI_Scatter(arr, 2, MPI_LONG_LONG,
 //            	rank_arr, 2, MPI_LONG_LONG, 0,
 //            	MPI_COMM_WORLD);
	// don't use this coz will use up a single rank's memory


	//int allocation = 1073741824 / my_mpi_size;
	int allocation = 1024 * 1024 / my_mpi_size;
	unsigned long long * rank_arr = calloc(allocation, sizeof(unsigned long long));
	unsigned long long rank_i = 0;

	// allocate arrays to individual ranks directly
	for(int i=0; i < allocation; i ++){
		rank_arr[i] = (allocation * my_mpi_rank) + (unsigned long long) i;
	}
	for(int i = 0; i < allocation; i ++){
		rank_i += rank_arr[i];
	}

	MPI_Barrier(MPI_COMM_WORLD); 

	start_cycles= GetTimeBase(); // ONLY EFFECTIVE AT BG/Q MACHINE
	/* reduce is called from every rank*/ 
	MPI_P2P_reduce(&rank_i, &sum, 1,
                       MPI_LONG_LONG, MPI_SUM, 0,
                       MPI_COMM_WORLD); 
	end_cycles= GetTimeBase();
	time_in_secs = ((double)(end_cycles - start_cycles)) / processor_frequency;

	MPI_Barrier(MPI_COMM_WORLD); 

	/* compare with MPI_reduce */
	start_cycles= GetTimeBase();
	MPI_Reduce(&rank_i, &sum_mpi, 1,
                   MPI_LONG_LONG, MPI_SUM, 0,
                   MPI_COMM_WORLD); 
	end_cycles= GetTimeBase();
	time_in_secs_mpi = ((double)(end_cycles - start_cycles)) / processor_frequency;
	/**/

	if(my_mpi_rank == 0){
		printf("%lld %.20f\n%lld %.20f\n", 
	 		sum, time_in_secs, sum_mpi, time_in_secs_mpi);
	}

	MPI_Barrier(MPI_COMM_WORLD); 

	#ifdef DEBUG
	if(my_mpi_rank == 0){
		printf("rank 0: sum %lld\n", sum);
		/* sum up the * sum */

	}
	#endif

	MPI_Finalize();

	free(rank_arr);


}

