/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/*   mam6@rpi.edu            **(*****************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "clcg4.h"

#include <mpi.h>
#include <pthread.h>
#include <unistd.h>


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

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

int N = 4;
int thread_n  = 4;
int ticks = 256;

#ifdef DEBUG
int hanging_tids = 0; //DEBUG
#endif

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_universeUpdate = PTHREAD_MUTEX_INITIALIZER;
#ifdef DEBUG
pthread_mutex_t lock_DEBUG = PTHREAD_MUTEX_INITIALIZER;
#endif

pthread_barrier_t barrier; 
pthread_barrier_t barrier_1; 
// pthread_barrierattr_t attr;

MPI_Request request_send_1;
MPI_Request request_send_2;
MPI_Request request_recv_1;
MPI_Request request_recv_2;

typedef struct update_args
{
    int ticks;
    int thread_no;
    int no_of_threads;
    int rows;
    int * sub_universe; 
    int * ALIVE_cells;
}update_arg;


void * update( void * args_ );


int main(int argc, char *argv[])
{
//    int i = 0;
    int mpi_myrank;
    int mpi_commsize;
    unsigned long long start_cycles=0; 
    unsigned long long end_cycles=0; 
    int time_in_secs_mpi;


// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 32,768 RNG streams - each row  has an independent *stream* -> and not a single value
/* value return by the RNG is greater than a set THRESHOLD, 
    then perform the above described basic rules for each cell. Otherwise, randomly pick state ofLIVEorDEAD*/
    
    InitDefault();

// Note, used the i to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
//     printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
//        mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    

// Allocate My rank’s chunk of the universe + space for "ghost" rows.
    int rank_rows = (int)(N / mpi_commsize);

    // int ** sub_universe = calloc(rank_rows + 2, sizeof(int *));
    // for(int i = 0; i < rank_rows + 2; i ++){
    //     sub_universe[i] = calloc(N, sizeof(int));
    // }

    // int sub_universe[rank_rows + 2][N];
    int * sub_universe = calloc((rank_rows + 2) * N, sizeof(int)); // a single row
    /*error: variable-sized object may not be initialized */
    /* use this trick: https://stackoverflow.com/questions/3082914/c-compile-error-variable-sized-object-may-not-be-initialized*/

    // memset( sub_universe, 0, (rank_rows + 2) * N * sizeof(int)); /* row-majored; sequential memory allocation; for mpi scatter*/
    // still supports indexing

    // int ** universe = calloc(N, sizeof(int *));
    // for(int i = 0; i < N; i ++){
    //     universe[i] = calloc(N, sizeof(int));
    // int universe[N][N];
    int * universe = calloc(N * N, sizeof(int));  // a single row
    // memset( universe, 0, N * N * sizeof(int)); // still supports indexing       

    start_cycles = GetTimeBase();


    // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //         void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    //         MPI_Comm comm)

    MPI_Scatter(universe, rank_rows * N, MPI_INT, sub_universe + N,  // ghost rows
                                rank_rows * N, MPI_INT, 0, MPI_COMM_WORLD);

    free(universe); 
    

// Initalize the universe (including "ghost" rows) with every cell being ALIVE
    for(int i = 0; i < rank_rows + 2; i ++){
        for(int j = 0; j < N; j ++){
            *(sub_universe + (i * N + j)) = ALIVE;
        }
    }    

// Create Pthreads here.
    pthread_t tid[thread_n];
    pthread_barrier_init(&barrier, NULL, thread_n);
    pthread_barrier_init(&barrier_1, NULL, thread_n);

    int * ALIVE_cells = calloc(ticks, sizeof(int)); /* modified by threads */
    int * ALIVE_cells_sum = calloc(ticks, sizeof(int));


    update_arg * args_ = calloc(thread_n, sizeof(update_arg));
    
    /* typedef struct update_args{
            int ticks;
            int thread_no;
            int no_of_threads;
            int rows; // no. of rows per rank (include 2 ghost rows)
            int ** sub_universe; 
        }update_arg;
    */
    #ifdef DEBUG
    hanging_tids = thread_n;
    #endif


    int tmp;

    for(int i = 0; i < thread_n; i ++){
        args_[i].ticks = ticks;
        args_[i].thread_no = i;
        args_[i].no_of_threads = thread_n;
        args_[i].rows = rank_rows + 2;
        args_[i].sub_universe = sub_universe; // thread share the rank's sub_universe and operate on it in parallel
        args_[i].ALIVE_cells = ALIVE_cells; // thread modify time ticks ALIVE_cells using mutex lock
        
        pthread_create(&tid[i], NULL, update, (void *)&args_[i]);

        /* int pthread_create(pthread_t *thread, const pthread_attr_t *attr,
                          void *(*start_routine) (void *), void *arg); 
            arg is passed as the sole argument for start_rountine;
                          */
    }

    for(int i = 0; i < thread_n; i ++){
        pthread_join(tid[i], (void *)&tmp);
    }

    
    

    // for(int i = 0; i < thread_n; i ++){

    //     #ifdef DEBUG
    //     while(1){
    //         printf("%d threads hanging. \n", hanging_tids);
    //         sleep(1);
    //     }
    //     #endif

    //     pthread_join(tid[i], (void *)&tmp); //todo; join seems not working
        
    //     #ifdef DEBUG
    //     printf("thread %d joined \n", i);
    //     #endif
    // }

    /* int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                      MPI_Datatype datatype, MPI_Op op, int root,
                      MPI_Comm comm) */
    MPI_Reduce(ALIVE_cells, ALIVE_cells_sum, ticks, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


    if(mpi_myrank == 0){
        end_cycles = GetTimeBase();
        time_in_secs_mpi = ((double)(end_cycles - start_cycles)) / g_processor_frequency;

    }
    /* some experiment required*/

    printf("time: %d, \n alive: ", time_in_secs_mpi);
    for(int i=0; i < ticks; i ++){
        printf("%d, ", ALIVE_cells[i]);
    }

    MPI_Barrier( MPI_COMM_WORLD );
    
    

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/



void * update( void * args_ ){

    /* within thread 
    
    args: 
    
    ticks
    thread no -> tid[i]:i
    no of threads per rank
    no of rows per *rank*  - include ghost rows
    sub_universe


    */

    update_arg * args = args_;

    int ticks = args->ticks;
    int thread_no = args->thread_no;
    int no_of_threads = args->no_of_threads;
    int rows = args->rows;
    int * sub_universe = args->sub_universe;  
    int * ALIVE_cells = args->ALIVE_cells;

    int mpi_myrank;
    int mpi_commsize;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);

    int alives;
    int * sub_universe_thread_part_copy = calloc(rows * N, sizeof(int)); /* rank_chunk * N */
    int living_nbrs;
    int state;

    /* notice thread using global time ticks*/

    /* 
        2. The expected heatmap for the 0% threshold is a completely dead universe except for maybe 
        a few cells alive at the right edge of the universe. If you don't see this, they you have some 
        sort of bug or other algorithm difference in your implementation. This happens because of the use 
        of a single copy universe which results in the initial 100% ALIVE universe turns to completely DEAD 
        after the first turn and never recovers due to the GOL rules.*/

    pthread_barrier_wait(&barrier_1);

    for( int t = 0; t < ticks; t++) { 
        alives = 0;

        /* Exchange row data with MPI ranks using MPI_Isend/Irecv from thread 0 
          w/i each MPI rank. 
          Yes, you must correctly MPI_Test or Wait to make sure messages operations correctly complete.
        
        [Note: have only 1 MPI rank/pthread perform ALL MPI operations per rank/thread group. 
        Dont’ allow multiple threads to perform any MPI operations within MPI rank/thread group.

        */
        if(thread_no == 0){ /* thread 0 perform MPI send/recv */

            #ifdef DEBUG
            printf("thread %d: here, start MPI send/recv\n", thread_no);
            #endif
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
                MPI_Isend(sub_universe + (1 * N), N, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &request_send_1);

            }else{
                MPI_Isend(sub_universe + (1 * N), N, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &request_send_1);   
            }
            
            /* int MPI_Wait(MPI_Request *request, MPI_Status *status)*/
            
            // MPI_Wait(&request_send_1, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated
            /* another way is to do while loop for MPI_Test and check `flag` */

            // #ifdef DEBUG
            // printf("thread %d: here, 0 \n", thread_no); // not reached
            // #endif


            if(mpi_myrank == mpi_commsize-1){
                MPI_Isend(sub_universe + (rows-2) * N, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &request_send_2);
            }else{
                MPI_Isend(sub_universe + (rows-2) * N, N, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &request_send_2);

            }

            // MPI_Wait(&request_send_1, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated

            // #ifdef DEBUG
            // printf("thread %d: here, 1 \n", thread_no); // reached
            // #endif
            if(mpi_myrank == 0){
                MPI_Irecv(sub_universe, N, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &request_recv_1);
            }else{
                MPI_Irecv(sub_universe, N, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &request_recv_1);
            }


            // MPI_Wait(&request_recv_1, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated

            // #ifdef DEBUG
            // printf("thread %d: here, 2 \n", thread_no); // reached
            // #endif
            if(mpi_myrank == mpi_commsize-1){
                MPI_Irecv(sub_universe + (rows - 1) * N, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &request_recv_2);
            }else{
                MPI_Irecv(sub_universe + (rows - 1) * N, N, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &request_recv_2);
            }

            // MPI_Wait(&request_recv_2, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated

            MPI_Wait(&request_send_1, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated
                #ifdef DEBUG
                printf("thread %d: here, 0 \n", thread_no); // not reached
                #endif
                        
            MPI_Wait(&request_send_2, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated
                #ifdef DEBUG
                printf("thread %d: here, 1 \n", thread_no); // reached
                #endif
            
            MPI_Wait(&request_recv_1, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated
                #ifdef DEBUG
                printf("thread %d: here, 2 \n", thread_no); // reached
                #endif
            
            MPI_Wait(&request_recv_2, MPI_STATUS_IGNORE); // Wait block until request succeed, `request` deallocated
                #ifdef DEBUG
                printf("thread %d: here, 3 \n", thread_no); // reached
                #endif

            // #ifdef DEBUG
            // printf("thread %d: here, 3 \n", thread_no); // reached
            // #endif
        }

        #ifdef DEBUG
        printf("thread %d: wait... \n", thread_no);
        #endif
        pthread_barrier_wait(&barrier); // thread level synchronization // todo: stuck here
        MPI_Barrier(MPI_COMM_WORLD); // rank level synchronization

        // #ifdef DEBUG
        // printf("rank %d: finished initialization. \n", mpi_myrank);
        // #endif


        /* HERE each PTHREAD can process a (chunk of?) * row *:

         - update universe making sure to use the correct row RNG stream 

         - factor in Threshold percentage as described 

         - use the right "ghost" row data at rank boundaries 

         - keep track of total number of ALIVE cells per tick across all threads w/i a MPI rank group. 

         
        */

        /* update *simutanously* */

        #ifdef DEBUG
        printf("thread %d: here, 0\n", thread_no);

        #endif

        /* get global(in rank) thread row index 0*/
        int thread_chunk = (int)((rows - 2) / no_of_threads); // exclude the ghost rows
        int thread_rows_0 = thread_chunk * (no_of_threads - 1) + 1; // include the ghost row



        /**RULES:

        Any live cell with fewer than two live neighbors dies, as if caused by under-population.
        Any live cell with two or three live neighbors lives on to the next generation.
        Any live cell with more than three live neighbors dies, as if by over-population.
        Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction.

        **/

        int get_state(int previous_state, int living_nbrs, int local_row, int rank_no, int rank_total){
            /*
            Each row has its own RNG stream
            */

            int global_row = rank_no * (N / rank_total) + local_row;
            long double rand_clcg = GenVal(global_row);
            double rand_rand = (double)rand() / (double)RAND_MAX; // [0, 1aq]

            if(rand_clcg > 0.25){ //
                
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

            }else{
                if(rand_rand > 0.5){
                    return ALIVE;
                }else{
                    return DEAD;
                }

            }

            //todo, use RNG

            return -1;

        }


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
            for(int j = 0; j < N; j ++){
                living_nbrs = 0;

                /* 2: previous row corresponding cell */
                living_nbrs += *(sub_universe + (thread_rows_0 + i - 1) * N + j);
                /* 6: next row corresponding cell */
                living_nbrs += *(sub_universe + (thread_rows_0 + i + 1) * N + j);

                if((j == 0) || (j == N-1)){ /* row ends; utilizing the end of the rows for missing left/right nbrs */
                    if(j == 0){
                        /* 1: previous row last cell*/ 
                        living_nbrs += *(sub_universe + (thread_rows_0 + i - 1) * N + (N-1)); //
                        /* 3: previous next cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i - 1) * N + (j+1)); // todo
                        /* 4: current row next cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i) * N + (j+1));
                        /* 5: next row next cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i + 1) * N + (j+1));
                        /* 7: next row last cell */ 
                        living_nbrs += *(sub_universe + (thread_rows_0 + i + 1) + (N-1)); //
                        /* 8: current row last cell */ 
                        living_nbrs += *(sub_universe + (thread_rows_0 + i) * N + (N+1)); //
                    }
                    else{
                        /* 1: previous row previous cell*/
                        living_nbrs += *(sub_universe + (thread_rows_0 + i - 1) * N + (j-1));
                        /* 3: previous head cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i - 1) * N); //
                        /* 4: current row head cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i) * N); //
                        /* 5: next row head cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i + 1) * N); //
                        /* 7: next row previous cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i + 1) * N + j-1);
                        /* 8: current row previous cell */
                        living_nbrs += *(sub_universe + (thread_rows_0 + i) * N + j-1);

                    }

                }else{
                    /* 1: previous row previous cell*/
                    living_nbrs += *(sub_universe + (thread_rows_0 + i - 1) * N + j-1);
                    /* 3: previous next cell */
                    living_nbrs += *(sub_universe + (thread_rows_0 + i - 1) * N + j+1);
                    /* 4: current row next cell */
                    living_nbrs += *(sub_universe + (thread_rows_0 + i) * N + j+1);
                    /* 5: next row next cell */
                    living_nbrs += *(sub_universe + (thread_rows_0 + i + 1) * N + j+1);
                    /* 7: next row previous cell */
                    living_nbrs += *(sub_universe + (thread_rows_0 + i + 1) * N + j-1);
                    /* 8: current row previous cell */
                    living_nbrs += *(sub_universe + (thread_rows_0 + i) * N + j-1);
                }

            // int get_state(int previous_state, int living_nbrs, int local_row, int rank_no, int rank_total){
                state = get_state(*(sub_universe + (thread_rows_0 + i) * N + j), living_nbrs, i, mpi_myrank, mpi_commsize);
                /* thread_rows index from non-ghost row*/
                *(sub_universe_thread_part_copy + (thread_rows_0 + i) * N + j) = state;

                if(state == ALIVE){
                    alives += 1;
                }


            }

        }

        #ifdef DEBUG
        printf("thread %d: here, 1\n", thread_no);

        #endif


        /*
        - keep track of total number of ALIVE cells per tick across all threads w/i a MPI rank group.

        - use pthread_mutex_trylock around shared counter variables **if needed**.

        */ 

        pthread_mutex_lock( &lock);
        ALIVE_cells[t] += alives; /*i : index of ticks*/
        pthread_mutex_unlock( &lock);

        // change the rank's sub_universe (rows * N) (in parallel)

        pthread_mutex_lock( &lock_universeUpdate );
        for(int i = 0; i < thread_chunk; i ++ ){ /* thread local row index*/
            for(int j = 0; j < N; j ++){
                *(sub_universe + (thread_rows_0 + i + 1) * N + j) = *(sub_universe_thread_part_copy + (thread_rows_0 + i + 1) * N + j); 
            }
        }
        pthread_mutex_unlock( &lock_universeUpdate );


    } /*tick end*/




    pthread_detach(pthread_self()); // nothing to return
    #ifdef DEBUG
    pthread_mutex_lock(&lock_DEBUG);
    hanging_tids -=1;
    pthread_mutex_unlock(&lock_DEBUG);
    #endif
    int * tmp = calloc(1, sizeof(int));
    *tmp = -1;
    return tmp;

}

