/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/*   mam6@rpi.edu            **(*****************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include<clcg4.h>

#include<mpi.h>
#include<pthread.h>

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime            
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

// You define these

int N = 1024;
int thread_n  = 64;
int ticks = 256;


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these


/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
    int mpi_myrank;
    int mpi_commsize;


// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 32,768 RNG streams - each * rank *  has an independent stream

    InitDefault();

// Note, used the i to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
       mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    

// Allocate My rank’s chunk of the universe + space for "ghost" rows.
    int ** sub_universe = calloc((int)(N / mpi_commsize) + 2, sizeof(int *));
    for(int i = 0; i < (int)(N / mpi_commsize) + 2; i ++){
        sub_universe[i] = calloc((int)(N / mpi_commsize), sizeof(int));
    }
    if(mpi_myrank == 0){
        int ** universe = calloc(N, sizeof(int *));
        for(int i = 0; i < N; i ++){
            universe[i] = calloc(N, sizeof(int));

         GetTimeBase();
    }
    

    // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //         void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    //         MPI_Comm comm)

        MPI_Scatter(universe, (int) (N / mpi_commsize), MPI_INT *, sub_universe+1,  // ghost rows
                                (int) (N / mpi_commsize), MPI_INT *, 0, MPI_COMM_WORLD)
    }

// Initalize the universe (including "ghost" rows) with every cell being ALIVE
    for(int i = 0; i < (int)(N / mpi_commsize) + 2; i ++){
        for(int j = 0; j < (int)(N / mpi_commsize); j ++){
            sub_universe[i][j] = ALIVE;
        }
    }    

// Create Pthreads here. All threads should go into the for-loop.
    pthread_t tid[thread_n];

    update_arg * args_ = calloc(thread_n, sizeof(update_arg));
    
    for(int i = 0; i < thread_n; i ++){
        args_[i].ticks = ;
        /*todo*/

        pthread_create(&tid[i], NULL, update, (void *) &args_[i]);
    }

    for(int i = 0; i < thread_n; i ++){
        pthread_join(tid[i], ALIVE_cells);
    }

    MPI_Reduce(ALIVE_cells_sujm, ALIVE_cells); /* todo*/





    MPI_Barrier( MPI_COMM_WORLD );
    
    

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/


void * update( void ** args_ ){

    /* within thread 
    
    args: 
    
    ticks
    thread no -> local
    no of threads per rank
    no of rows per *rank*
    sub_universe


    */

    update_arg * args = args_;

    int ticks = args->ticks;
    int thread_no = args->thread_no;
    int no_of_threads = args->no_of_threads;
    int rows = args->rows;
    int ** sub_universe = args->sub_universe;  

    int mpi_myrank;
    int mpi_commsize;
    MPI_Request request;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);

    for( i = 0; i < ticks; i++) { 

        /* Exchange row data with MPI ranks using MPI_Isend/Irecv from thread 0 
          w/i each MPI rank. 
          Yes, you must correctly MPI_Test or Wait to make sure messages operations correctly complete.
        
        [Note: have only 1 MPI rank/pthread perform ALL MPI operations per rank/thread group. 
        Dont’ allow multiple threads to perform any MPI operations within MPI rank/thread group.

        */
        if(thread_no == 0){
            /* rank 0 <-> rank 1 

                rank 0: row -1     <- rank (total - 1): row END
                        row END+1  <- rank 1: row 0

                rank 1: row -1     <- rank 0: row END
                        row END+1  <- rank 2: row 0

                ...

                rank (total -1 ): row -1     <- rank (total - 2): row END
                                  row END+1  <- rank 0: row 0

                each rank: 2 * recv, 2 * send

            */

            
            /* int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
                int tag, MPI_Comm comm, MPI_Request *request) */
            /* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
               int source, int tag, MPI_Comm comm, MPI_Request *request) */
            if(mpi_myrank == 0){
                MPI_Isend(sub_universe[1], N, MPI_INT, mpi_commsize-1, MPI_tag, MPI_COMM_WORLD, &request);
            }else{
                MPI_Isend(sub_universe[1], N, MPI_INT, mpi_myrank - 1, MPI_tag, MPI_COMM_WORLD, &request);   
            }
            if(mpi_myrank == mpi_commsize-1){
                MPI_Isend(sub_universe[rows-1], N, MPI_INT, 0, MPI_tag, MPI_COMM_WORLD, &request);
            }else{
                MPI_Isend(sub_universe[rows-1], N, MPI_INT, mpi_myrank + 1, MPI_tag, MPI_COMM_WORLD, &request);

            }

            while(1){
                MPI_Irecv();
                if() break; /*todo*/
            }
            while(1){
                MPI_Irecv();
                if() break; /*todo*/
            }

        }


        /* HERE each PTHREAD can process a (chunk of?) * row *:

         - update universe making sure to use the correct row RNG stream 

         - factor in Threshold percentage as described 

         - use the right "ghost" row data at rank boundaries 

         - keep track of total number of ALIVE cells per tick across all threads w/i a MPI rank group. 

         
        */

        /* update *simutanously* */
        /* mind the corners ! -- asked on submitty */
        int ** sub_universe_copy;
        int thread_chunk = (int)(rows / thread_no);

        /*todo*/
        /*****/


        /*
        - keep track of total number of ALIVE cells per tick across all threads w/i a MPI rank group.

        - use pthread_mutex_trylock around shared counter variables **if needed**.

        */ 

        ALIVE_cells[i] += alives;

        /*todo*/


    } /*tick end*/


    return &ALIVE_cells; // each rank; int[ticks] -> no. of alive cells at each time tick

}


