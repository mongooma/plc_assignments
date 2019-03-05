#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
/* */
#include "functions.h"


// int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
//                       MPI_Datatype datatype, MPI_Op op, int root,
//                       MPI_Comm comm)

/* A Code snippet:

	#define  LEN   1000

	float val[LEN];        // local array of values 
	int count;             // local number of values 
	int myrank, minrank, minindex;
	float minval;

	struct {
	   float value;
	   int   index;
	} in, out;

	//local minloc 
	in.value = val[0];
	in.index = 0;
	for (i=1; i < count; i++)
	   if (in.value > val[i]) {
	       in.value = val[i];
	       in.index = i;
	   }

	// global minloc 

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	in.index = myrank*LEN + in.index;
	MPI_Reduce( in, out, 1, MPI_FLOAT_INT, MPI_MINLOC, root, comm );
	// At this point, the answer resides on process root
	    
	if (myrank == root) {
	// read answer out
	    
	   minval = out.value;
	   minrank = out.index / LEN;
	   minindex = out.index % LEN;

*/


int MPI_P2P_reduce(const int *sendbuf, int * recvbuf, int count,
                     MPI_Datatype datatype, MPI_Op op, int root,
                     MPI_Comm comm){

	/*  your implementation will only perform the MPI SUM operation 
		and the ﬁnal reduction result goes to MPI rank 0. */
	/* 	Take in the same params as MPI_reduce 
	 	
	 	The  global  reduce  functions  (MPI_Reduce,  MPI_Op_create,  MPI_Op_free,  MPI_Allreduce, MPI_Reduce_scatter,
       	MPI_Scan) perform a global reduce operation (such as sum, max, logical AND, etc.) across all the members of  a
       	group.  The reduction operation can be either one of a predefined list of operations, or a user-defined opera‐
       	tion.*/

	/* 	params:
		
		(int)
		sendbuf   Address of send buffer (choice).

	    count     Number of elements in send buffer (integer).

	    datatype  Data type of elements of send buffer (handle).

	    op        Reduce operation (handle). //
	    The reduction functions ( MPI_Op ) do not return an error value

	    The  following  predefined  operations  are  supplied  for  MPI_Reduce  and  related  functions MPI_Allreduce,
       	MPI_Reduce_scatter, and MPI_Scan. These operations are invoked by placing the following in op:

            Name                Meaning
            ---------           --------------------
            MPI_MAX             maximum
            MPI_MIN             minimum
            MPI_SUM             sum
            MPI_PROD            product
            MPI_LAND            logical and
            MPI_BAND            bit-wise and
            MPI_LOR             logical or
            MPI_BOR             bit-wise or
            MPI_LXOR            logical xor
            MPI_BXOR            bit-wise xor
            MPI_MAXLOC          max value and location
            MPI_MINLOC          min value and location

			 Op                       Allowed Types
            ----------------         ---------------------------
            MPI_MAX, MPI_MIN         C integer, Fortran integer,
                                     floating-point

            MPI_SUM, MPI_PROD        C integer, Fortran integer,
                                     floating-point, complex

            MPI_LAND, MPI_LOR,       C integer, logical
            MPI_LXOR

            MPI_BAND, MPI_BOR,       C integer, Fortran integer, byte
            MPI_BXOR


 		For MPI_MAXLOC / MINLOC: (return the value and the index of the rank )

           Name                 Description
           MPI_FLOAT_INT            float and int
           MPI_DOUBLE_INT           double and int
           MPI_LONG_INT             long and int
           MPI_2INT                 pair of ints
           MPI_SHORT_INT            short and int
           MPI_LONG_DOUBLE_INT      long double and int

	    root      Rank of root process (integer).

	    comm      Communicator (handle).
	    
	    (output)
	    recvbuf   Address of receive buffer (choice, significant only at root)..
		
		Returns:
		

	*/


	/*1.  Each rank computes sum over local data array.*/

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	for(int i = 0; i < count; i ++){
		*recvbuf += sendbuf[i];
	}
	#ifdef DEBUG
		printf("rank %d: local sum %d. \n ", myrank, *recvbuf);
	#endif


	/*2. Recursively perform pairwise sums for a higher rank to a lower rank 
		Using MPI_Isend/Irecv 
	*/

	int r = 0;
	int step = 1;
	int rev = 0;
	int *tmp_buf = calloc(1, sizeof(int));

	MPI_Request request;

	MPI_Barrier(MPI_COMM_WORLD); /* synchronization */
	while(1){
		// step = (int) pow(2, (double) r);
		rev = 0;
		/* from high to low*/
		for(int i = atoi(getenv("RANK_NO")) - step; i >= 0; i -= step){
			/* This is supposed to be very slow because we are using for loop here*/
			if(!rev){
				/*
				int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
            	int tag, MPI_Comm comm, MPI_Request *request)
            	*/
            	if(myrank == i){
					MPI_Isend(recvbuf, 1, MPI_INT, i - step, 0, MPI_COMM_WORLD, &request);
				#ifdef DEBUG
					printf("rank %d: send %d to rank %d. \n ", myrank, *recvbuf, i - step);
				#endif
            	}

			}else{
				if(myrank == i){
					MPI_Irecv(tmp_buf, 1, MPI_INT, i + step, 0, MPI_COMM_WORLD, &request);
				#ifdef DEBUG
					printf("rank %d: recv %d from rank %d. \n", myrank, *tmp_buf, i + step);
				#endif
					/* sum up */
					*recvbuf += *tmp_buf;

				}
				
			}

			MPI_Barrier(MPI_COMM_WORLD); /* use barrier to synchronize */

			rev = !rev;
		}

		r++;
		step = step * 2;
		// if((int) pow( 2, r) > atoi(getenv("RANK_NO"))) break; // All the results in rank 0
		if((step) >= atoi(getenv("RANK_NO"))) break; // All the results in rank 0
	}

	/* The result will be stored at rank 0 *sum*/
	free(tmp_buf);
	return 0;
}

