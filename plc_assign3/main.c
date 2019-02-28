#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BGQ 1 // when running BG/Q, comment out when running on mastiff
#ifdef BGQ 
#include<hwi/include/bqc/A2_inlines.h> 
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


int main(){

	/* Use the Blue Gene/Q’s “GetTimeBase()” function to measure the number of 
		cycles this point-2-point reduce operation took as well as the MPI 
		collection MPI reduce
	*/
	double time_in_secs = 0; 
	double processor_frequency = 1600000000.0; 
	unsigned long long start_cycles=0; 
	unsigned long long end_cycles=0; 

	/*...*/

	start_cycles= GetTimeBase(); 
	MPI_P2P_Reduce(); 
	end_cycles= GetTimeBase();

	time_in_secs = ((double)(end_cycles - start_cycles)) / processor_frequency;




	printf("%d %d\n%d %d\n", 
			sum_user, time_in_secs_user, sum, time_in_secs)


}

