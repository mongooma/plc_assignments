pthread_barrier_t barrier; 
pthread_barrierattr_t attr;

MPI_Request request_send_1;
MPI_Request request_send_2;
MPI_Request request_recv_1;
MPI_Request request_recv_2;

#ifdef DEBUG
int hanging_tids = 0; //DEBUG
#endif

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_universeUpdate = PTHREAD_MUTEX_INITIALIZER;
#ifdef DEBUG
pthread_mutex_t lock_DEBUG = PTHREAD_MUTEX_INITIALIZER;
#endif


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


