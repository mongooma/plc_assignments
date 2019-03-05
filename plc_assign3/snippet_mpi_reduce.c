/* code snippet for MPI_reduce */


/* Example 3: Each process has an array of 30 doubles, in C. For each of the 30 locations, compute  the  value
       and rank of the process containing the largest value. 
*/
/* each process has an array of 30 double: ain[30]
*/
double ain[30], aout[30];
int  ind[30];
struct {
   double val;
   int   rank;
} in[30], out[30];
int i, myrank, root;

MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
for (i=0; i<30; ++i) {
   in[i].val = ain[i];
   in[i].rank = myrank;
}
MPI_Reduce( in, out, 30, MPI_DOUBLE_INT, MPI_MAXLOC, root, comm );
/* At this point, the answer resides on process root
*/
if (myrank == root) {
   /* read ranks out
    */
   for (i=0; i<30; ++i) {
       aout[i] = out[i].val;
       ind[i] = out[i].rank;
   }
}