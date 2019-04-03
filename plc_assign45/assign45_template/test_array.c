#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include"clcg4.h"

// #include<mpi.h>
#include<pthread.h>

int main(){

	int l[2][2];
	memset(l, 0, 2 * 2 * sizeof(int));

	l[1][0] = 1; // still supports indexing

	printf("l[1][0]: %d \n", *(l + 1)[0]);

	return 0;
}