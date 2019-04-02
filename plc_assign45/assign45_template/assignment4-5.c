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
#include <hwi/include/bqc/A2_inlines.h>
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

typedef struct update_args
{
    int ticks;
    int thread_no;
    int no_of_threads;
    int rows;
    int ** sub_universe; 
}update_arg;

// You define these

int N = 1024;
int thread_n  = 64;
int ticks = 256;
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

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
    int rank_rows = (int)(N / mpi_commsize);
    int ** sub_universe = calloc(rank_rows + 2, sizeof(int *));
    for(int i = 0; i < rank_rows + 2; i ++){
        sub_universe[i] = calloc(N, sizeof(int));
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

        MPI_Scatter(universe, rank_rows, MPI_INT *, sub_universe+1,  // ghost rows
                                rank_rows, MPI_INT *, 0, MPI_COMM_WORLD)
    }

// Initalize the universe (including "ghost" rows) with every cell being ALIVE
    for(int i = 0; i < rank_rows + 2; i ++){
        for(int j = 0; j < N; j ++){
            sub_universe[i][j] = ALIVE;
        }
    }    

// Create Pthreads here.
    pthread_t tid[thread_n];
    int * ALIVE_cells = call(ticks, sizeof(int)); /* modified by threads */
    int * ALIVE_cells_sum = call(ticks, sizeof(int));


    update_arg * args_ = calloc(thread_n, sizeof(update_arg));
    
    /* typedef struct update_args{
            int ticks;
            int thread_no;
            int no_of_threads;
            int rows; // no. of rows per rank (include 2 ghost rows)
            int ** sub_universe; 
        }update_arg;
    */
    for(int i = 0; i < thread_n; i ++){
        args_[i].ticks = ticks;
        args_[i].thread_no = i;
        args_[i].no_of_threads = thread_n;
        args_[i].rows = rank_rows + 2;
        args_[i].sub_universe = sub_universe; // thread share the rank's sub_universe and operate on it in parallel
        args_[i].ALIVE_cells = ALIVE_cells; // thread modify time ticks ALIVE_cells using mutex lock
        
        pthread_create(&tid[i], NULL, update, (void *) &args_[i]);
    }

    for(int i = 0; i < thread_n; i ++){
        pthread_join(tid[i], NULL);
    }

    /* int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                      MPI_Datatype datatype, MPI_Op op, int root,
                      MPI_Comm comm) */
    MPI_Reduce(ALIVE_cells, ALIVE_cells_sum, ticks, MPI_INT, MPI_Sum, 0, MPI_COMM_WORLD);


    /* some experiment required*/


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
    thread no -> tid[i]:i
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
    int * ALIVE_cells = args->ALIVE_cells;

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
            /* int MPI_Irecv(      void *buf, int count, MPI_Datatype datatype, int source, 
                int tag, MPI_Comm comm, MPI_Request *request) */
            if(mpi_myrank == 0){
                MPI_Isend(sub_universe[1], N, MPI_INT, mpi_commsize - 1, MPI_tag, MPI_COMM_WORLD, &request);
            }else{
                MPI_Isend(sub_universe[1], N, MPI_INT, mpi_myrank - 1, MPI_tag, MPI_COMM_WORLD, &request);   
            }
            
            /* int MPI_Wait(MPI_Request *request, MPI_Status *status)*/
            MPI_Wait(&request, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated
            /* another way is to do while loop for MPI_Test and check `flag` */

            if(mpi_myrank == mpi_commsize-1){
                MPI_Isend(sub_universe[rows-1], N, MPI_INT, 0, MPI_tag, MPI_COMM_WORLD, &request);
            }else{
                MPI_Isend(sub_universe[rows-1], N, MPI_INT, mpi_myrank + 1, MPI_tag, MPI_COMM_WORLD, &request);

            }

            MPI_Wait(&request, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated

            if(mpi_myrank == 0){
                MPI_Irecv(sub_universe[0], N, MPI_INT, mpi_commsize - 1, MPI_tag, MPI_COMM_WORLD, &request);
            }else{
                MPI_Irecv(sub_universe[0], N, MPI_INT, mpi_myrank - 1, MPI_tag, MPI_COMM_WORLD, &request);
            }

            MPI_Wait(&request, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated

            if(mpi_myrank == mpi_commsize-1){
                MPI_Irecv(sub_universe[rows], N, MPI_INT, 0, MPI_tag, MPI_COMM_WORLD, &request);
            }else{
                MPI_Irecv(sub_universe[rows], N, MPI_INT, mpi_myrank + 1, MPI_tag, MPI_COMM_WORLD, &request);
            }

            MPI_Wait(&request, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated

        }


        /* HERE each PTHREAD can process a (chunk of?) * row *:

         - update universe making sure to use the correct row RNG stream 

         - factor in Threshold percentage as described 

         - use the right "ghost" row data at rank boundaries 

         - keep track of total number of ALIVE cells per tick across all threads w/i a MPI rank group. 

         
        */

        /* update *simutanously* */
        int ** sub_universe_thread_part_copy; /* thread_chunk * N */

        /* get global(in rank) thread row index 0*/
        int thread_chunk = (int)((rows - 2) / thread_no); // exclude the ghost rows
        int thread_rows_0 = thread_chunk * (thread_no - 1) + 1; // include the ghost row

        /**RULES:

        Any live cell with fewer than two live neighbors dies, as if caused by under-population.
        Any live cell with two or three live neighbors lives on to the next generation.
        Any live cell with more than three live neighbors dies, as if by over-population.
        Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction.

        **/

        int living_nbrs;
        int state;
        int get_state(int previous_state, int living_nbrs){

            if(living_nbrs < 2){
                return DEAD;
            }else if (living_nbrs == 2)
            {
                if(previous_state == ALIVE){
                    return ALIVE;
                }else{
                    return DEAD;
                }
            }else if(living_nbrs == 3){
                return ALIVE;
            }else if(living_nbrs > 3){
                return DEAD;
            }

        }

        int alives;

        for(int i = 0; i < thread_chunk; i ++ ){ /* thread local row index*/
            /*
            +++++++++++++++++++++++++++++++ row = 0 (rank ghost row)
            ...............................
                      ...
            ...............................
            ...............................  i = 1
                      ...
            ...............................  i = thread_chunk
            ...............................
                      ...
            ...............................  
            +++++++++++++++++++++++++++++++ row = rows (rank ghost row)


            123
            8_4
            765  ---> eight nbrs

            */
            for(int j = 0: j < N; j ++){
                living_nbrs = 0

                /* 2: previous row corresponding cell */
                living_nbrs += sub_universe[thread_rows_0 + i - 1][j];
                /* 6: next row corresponding cell */
                living_nbrs += sub_universe[thread_rows_0 + i + 1][j];

                if((j == 0) || (j == N-1)){ /* row ends; utilizing the end of the rows for missing left/right nbrs */
                    if(j == 0){
                        /* 1: previous row last cell*/ 
                        living_nbrs += sub_universe[thread_rows_0 + i - 1][N-1]; //
                        /* 3: previous next cell */
                        living_nbrs += sub_universe[thread_rows_0 + i - 1][j+1];
                        /* 4: current row next cell */
                        living_nbrs += sub_universe[thread_rows_0 + i ][j+1];
                        /* 5: next row next cell */
                        living_nbrs += sub_universe[thread_rows_0 + i + 1][j+1];
                        /* 7: next row last cell */ 
                        living_nbrs += sub_universe[thread_rows_0 + i + 1][N-1]; //
                        /* 8: current row last cell */ 
                        living_nbrs += sub_universe[thread_rows_0 + i][N+1]; //
                    }
                    else{
                        /* 1: previous row previous cell*/
                        living_nbrs += sub_universe[thread_rows_0 + i - 1][j-1];
                        /* 3: previous head cell */
                        living_nbrs += sub_universe[thread_rows_0 + i - 1][0]; //
                        /* 4: current row head cell */
                        living_nbrs += sub_universe[thread_rows_0 + i ][0]; //
                        /* 5: next row head cell */
                        living_nbrs += sub_universe[thread_rows_0 + i + 1][0]; //
                        /* 7: next row previous cell */
                        living_nbrs += sub_universe[thread_rows_0 + i + 1][j-1];
                        /* 8: current row previous cell */
                        living_nbrs += sub_universe[thread_rows_0 + i][j-1];

                    }

                }else{
                    /* 1: previous row previous cell*/
                    living_nbrs += sub_universe[thread_rows_0 + i - 1][j-1];
                    /* 3: previous next cell */
                    living_nbrs += sub_universe[thread_rows_0 + i - 1][j+1];
                    /* 4: current row next cell */
                    living_nbrs += sub_universe[thread_rows_0 + i ][j+1];
                    /* 5: next row next cell */
                    living_nbrs += sub_universe[thread_rows_0 + i + 1][j+1];
                    /* 7: next row previous cell */
                    living_nbrs += sub_universe[thread_rows_0 + i + 1][j-1];
                    /* 8: current row previous cell */
                    living_nbrs += sub_universe[thread_rows_0 + i][j-1];
                }

                state = get_state(sub_universe[thread_rows[0] + i + 1][j], living_nbrs);
                /* thread_rows index from non-ghost row*/
                sub_universe_thread_part_copy[thread_rows[0] + i + 1][j] = state;

                if(state == ALIVE){
                    alives += 1;
                }


            }

        }


        /*
        - keep track of total number of ALIVE cells per tick across all threads w/i a MPI rank group.

        - use pthread_mutex_trylock around shared counter variables **if needed**.

        */ 

        pthread_mutex_lock( &lock);
        ALIVE_cells[i] += alives; /*i : index of ticks*/
        pthread_mutex_unlock( &lock);


    } /*tick end*/


    return NULL;

}


