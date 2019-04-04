

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


            // #ifdef DEBUG
            // printf("thread %d: here, 3 \n", thread_no); // reached
            // #endif

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
        }

        pthread_barrier_wait(&barrier); // thread level synchronization
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

            if(rand_clcg > 0.25){
                
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