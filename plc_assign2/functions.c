#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>
#include "functions.h"

/* this file will be accessed per MPI rank, 
   the bins used by each rank are stored in bin_rank
*/

#define block_size 32 /* not a variable in this assignment */
#define c_minus_1 0 /*primary value to set off the calculation chain*/

char lookup[16][5] = { "0000",
						"0001",
						"0010",
						"0011",
						"0100",
						"0101",
						"0110",
						"0111",
						"1000",
						"1001",
						"1010",
						"1011",
						"1100",
						"1101",
						"1110",
						"1111" };


/*notation: significance order: b<--|<--a*/




void readInData(const char * filename,  
						char * hex_input_a, char * hex_input_b){
	/* to debug */

	FILE *my_input_file=NULL;

	#ifdef DEBUG
	#endif

	if( (my_input_file = fopen( filename, "r")) == NULL ){
		printf("Failed to open input data file: %s \n", filename);
		exit(-1);
	}

	fscanf( my_input_file, "%s %s", hex_input_a, hex_input_b );

	#ifdef DEBUG
	#endif
	
	fclose( my_input_file );

}

int convert_hex_2_bit(char * hex_input_a, char * hex_input_b, int * bin1, int * bin2, int input_size){
	/* convert a whole hex array and store into two bins with whole length*/
	/* implement a look-up table instead */
	
	char *tmp1, *tmp2, tmp_[2];
	tmp_[1] = '\0';

	/* hex1 and hex2 of length input_size + 1*/
	for (int i = 0; i < input_size; i++) {
		tmp_[0] = hex_input_a[i];
		tmp1 = lookup[strtol(tmp_, 0, 16)];
		tmp_[0] = hex_input_b[i];
		tmp2 = lookup[strtol(tmp_, 0, 16)];
		for (int j = 0; j < 4; j++) {
			tmp_[0] = tmp1[j];
			bin1[i * 4 + j] = atoi(tmp_); /* atoi will return exact as it looks*/
			tmp_[0] = tmp2[j];
			bin2[i * 4 + j] = atoi(tmp_); /* atoi will return exact as it looks*/
		}
	}

	return 0;

}

int revert_binary(int * bin1, int * bin2, int bits){
	/* revert the binary string to set the last digit as the MSB*/
	for (int i = 0; i < (bits / 2); i++) {

		int tmp = bin1[i];
		bin1[i] = bin1[(bits - 1) - i];
		bin1[(bits - 1) - i] = tmp;
		tmp = bin2[i];
		bin2[i] = bin2[(bits - 1) - i];
		bin2[(bits - 1) - i] = tmp;
	}

	return 0; 

}


int revert_hex_sum(int * sumi, int bits){

	for (int i = 0; i < bits / 2; i++) {
		int tmp = sumi[i];
		sumi[i] = sumi[(bits - 1) - i];
		sumi[(bits - 1) - i] = tmp;
	}

	return 0; 

}

int convert_bit_2_hex(int * sumi, int bits){

	FILE *f = fopen("output.txt", "w"); /* use fprintf for debugging output*/
	for (int i = 0; i < bits; i += 4) {
		int tmp = 0; 
		for (int j = 0; j < 4; j++) {
			tmp += (sumi[i + j] % 2);
			if (j == 3) break;
			tmp = tmp * 2;
		}
		fprintf(f, "%X", tmp); /* %x format for hex*/
	}
	fprintf(f, "\n");

	fclose(f);

	return 0; 
}


/* Carry-Lookahead Adder */
/* master */
/* cla executed within one rank */
int cla(int use_barrier, int my_mpi_rank, int my_mpi_size, int alloc, int * bin1, int *bin2, int *sumi){
	/* bin1 and bin2 are local bin arrays for each rank, i.e. bin_rank_1, bin_rank_2*/	
	
	const int bits = alloc; /* every hex character is 4 digit binary */
	const int ngroups = bits/block_size; 
	const int nsections = ngroups/block_size; 
	const int nsupersections = nsections/block_size;


	/* local definitions of the various arrays used */
	int * gi = calloc(bits, sizeof(int)); 
	int * pi = calloc(bits, sizeof(int));  
	int * ci = calloc(bits, sizeof(int)); 

	int * ggj = calloc(ngroups, sizeof(int));
	int * gpj = calloc(ngroups, sizeof(int));
	int * gcj = calloc(ngroups, sizeof(int));

	int * sgk = calloc(nsections, sizeof(int));
	int * spk = calloc(nsections, sizeof(int));
	int * sck = calloc(nsections, sizeof(int));
	
	int * ssgl = calloc(nsupersections, sizeof(int));
	int * sspl = calloc(nsupersections, sizeof(int));
	int * sscl = calloc(nsupersections, sizeof(int));

	#ifdef DEBUG
		printf("cla: here0! \n"); /* three reached here */
	#endif 


	/* all ranks except Rank 0, should post an MPI Irecv message before the cla 
	calculations start */ 
	/* what inputs are they waiting for ?*/
	/* is this the c-minus for the left most group/section/nsection ?*/
	
	/* perform cla for each rank*/

	/* g */
	/* p */
	g(bits, gi, bin1, bin2);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	p(bits, pi, bin1, bin2);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}	

	/* between each algorithm step (level) perform an MPI Barrier collective operation to 
	keep all ranks in step with each other. */

	/* gg */
	/* gp */
	gg(ngroups, ggj, gi, pi);
	
	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}
	
	gp(ngroups, gpj, pi);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	/* sg */
	/* sp */
	sg(nsections, sgk, ggj, gpj);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	sp(nsections, spk, gpj);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	/* ssg */
	/* ssp */
	ssg(nsupersections, ssgl, sgk, spk);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	ssp(nsupersections, sspl, spk);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	/* ssc */
	ssc(sscl, ssgl, sspl, my_mpi_size, my_mpi_rank, nsupersections);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	/* sc */
	sc(sscl, sck, sgk, spk, nsections, my_mpi_rank);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	/* gc */

	gc(sck, gcj, ggj, gpj, ngroups, my_mpi_rank);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	/* c */

	c(gcj, ci, gi, pi, bits, my_mpi_rank);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	sum_cla(sumi, bin1, bin2, ci, bits);

	if(use_barrier){MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/}

	return 0;

}


/* Some lower-level function basics */


/*1st level - bits 4096 */
int g(int bits, int *gi, int *bin1, int * bin2){  /*generates; 1 1*/
	
	for (int i = 0; i < bits; i++) { 

		gi[i] = bin1[i] & bin2[i];
	}

	return 0;
}


int p(int bits, int *pi, int *bin1, int * bin2) { /*propagates: 1 0, 0 1, 1 1*/
	for (int i = 0; i < bits; i++) {

		pi[i] = bin1[i] | bin2[i];
	}
	return 0;
}

/*2nd level - group 512*/
int gg(int ngroups, int *ggj, int *gi, int * pi) {

	for (int j = 0; j < ngroups; j++) {
		for(int b = 0+1; b < block_size; b++){
			int tmp = 1;
			for(int b_minus = 0; b_minus < b; b_minus++){ 
				if (b_minus == b-1) { tmp &= gi[j * block_size + (block_size - 1 - b_minus)]; /* b_minus maximum is block_size - 1*/}
								else{ tmp &= pi[j * block_size + (block_size - 1 - b_minus)];}
			}
			ggj[j] |= tmp;
		}
	}
	return 0;

}

int gp(int ngroups, int *gpj, int * pi) { /*the logic will be the same for other propagate functions at higher levels:
			in order to propagate from previous group, the chain: 1-->2-->3-->4--> must be valid at every point
		   */
	for (int j = 0; j < ngroups; j++) {
		gpj[j] = 1;
		for (int b = block_size-1; b > -1; b--){
			gpj[j] &= pi[j * block_size + b];
		}
	}
	return 0;

}

/* 3rd level - sections 64*/

int sg(int nsections, int *sgk, int *ggj, int * gpj) { /*check gg()*/

	for (int k = 0; k < nsections; k++) {
		for(int b = 0+1; b < block_size; b++){
			int tmp = 1;
			for(int b_minus = 0; b_minus < b; b_minus++){ 
				if (b_minus == b-1) { tmp &= ggj[k * block_size + (block_size - 1 - b_minus)]; /* b_minus maximum is block_size - 1*/}
								else{ tmp &= gpj[k * block_size + (block_size - 1 - b_minus)];}
			}
			sgk[k] |= tmp;
		}
	}
	return 0;
}

int sp(int nsections, int *spk, int *gpj) {

	for (int k = 0; k < nsections; k++) {
		spk[k] = 1;
		for (int b = block_size-1; b > -1; b--){
			spk[k] &= gpj[k * block_size + b];
		}
	}
	return 0;

}

/* 4th level - supersection 8*/

int ssg(int nsupersections, int *ssgl, int *sgk, int * spk) { /*check gg()*/
	for (int l = 0; l < nsupersections; l++) {
		for(int b = 0+1; b < block_size; b++){
			int tmp = 1;
			for(int b_minus = 0; b_minus < b; b_minus++){ 
				if (b_minus == b-1) { tmp &= sgk[l * block_size + (block_size - 1 - b_minus)]; /* b_minus maximum is block_size - 1*/}
								else{ tmp &= spk[l * block_size + (block_size - 1 - b_minus)];}
			}
			ssgl[l] |= tmp;
		}
	}
	return 0;
}

int ssp(int nsupersections, int *sspl, int * spk) { 

	for (int l = 0; l < nsupersections; l++) {
		sspl[l] = 1;
		for (int b = block_size-1; b > -1; b--){
			sspl[l] &= spk[l * block_size + b];
		}
	}
	return 0;

}

/* reduce: carries computation
*/

int ssc(int * sscl, int * ssgl, int * sspl, int my_mpi_size, int my_mpi_rank, int nsupersections) { /*at supersection level*/
	/* use blocking msg passing for now*/

	/* here we communicate between ranks since we need sscl[l-1] from neighbor node*/
	/* Use MPI Isend and MPI Irecv routines to exchange the nearest-neighbor 
	carry bit information at the sscl level*/
	/* so that all processes share the information of ssc_{l-1}, l= k/block_size  
		for next step sc_{k-1} correction, k mod block_size = 0*/
	/* send and recv could be called globally since source is stated*/
	
	MPI_Request request;

	int buf;

	#ifdef DEBUG
		printf("rank %d, ssc: here0.\n", my_mpi_rank);
	#endif


	for(int i = 0; i < my_mpi_size; i++){

		if ((my_mpi_rank == i) && ((i == 0) || (my_mpi_size == 1)) ){

			buf = c_minus_1;
			for (int l = 0; l < nsupersections; l++) {

					if (l == 0) { sscl[l] = ssgl[l] | (sspl[l] & buf );}
							else{ sscl[l] = ssgl[l] | (sspl[l] & sscl[l-1]);}
			}

			if(my_mpi_size != 1){
				MPI_Isend(&(sscl[nsupersections-1]), 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &request); 
			}else{ return 0;}

		} 
		
		MPI_Barrier(MPI_COMM_WORLD); /*To have every process synchronized*/ 

		#ifdef DEBUG
			printf("rank %d, ssc: here1.\n", my_mpi_rank);
		#endif

		if( (my_mpi_rank == i) && ( i > 0) ){
			/* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
           int source, int tag, MPI_Comm comm, MPI_Request *request)*/
			MPI_Irecv(&buf, 1, MPI_INT, my_mpi_rank - 1, 0, MPI_COMM_WORLD, &request);

			for (int l = 0; l < nsupersections; l++) {

				if (l == 0) { sscl[l] = ssgl[l] | (sspl[l] & buf );}
						else{ sscl[l] = ssgl[l] | (sspl[l] & sscl[l-1]);}

			}

		}

		MPI_Barrier(MPI_COMM_WORLD); /*To have every process synchronized*/ 

		#ifdef DEBUG
			printf("rank %d, ssc: here2.\n", my_mpi_rank);
		#endif

		
		if((my_mpi_rank == i) && (i < my_mpi_size-1)){
			/* int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
			int tag, MPI_Comm comm, MPI_Request *request)*/
			MPI_Isend(&(sscl[nsupersections-1]), 1, MPI_INT, my_mpi_rank+1, 0, MPI_COMM_WORLD, &request); 
		}

		MPI_Barrier(MPI_COMM_WORLD); /*To have every process synchronized*/ 

		#ifdef DEBUG
			printf("rank %d, ssc: here3.\n", my_mpi_rank);
		#endif

	}

	return 0;

}

int sc(int *sscl, int * sck, int * sgk, int * spk, int nsections, int my_mpi_rank) { /*section level; for now independently computed at current level*/

	for (int k = 0; k < nsections; k++) {

		/* replacing each sck−1 when k mod 32 = 0 using correct s2cl−1,l = k div 32 */
		if ((k % block_size == 0) && (k != 0)){ /*the last one in the section*/
			sck[k-1] = sscl[k / block_size - 1]; /* integer division, i.e. 10/3 = 3*/
		}

		if ((my_mpi_rank == 0) && (k == 0)) {
			sck[k] = sgk[k] | (spk[k] & c_minus_1);
		}
		else {
			sck[k] = sgk[k] | (spk[k] & sck[k - 1]);
		}

	}

	return 0;
}

int gc(int * sck, int * gcj, int * ggj, int * gpj, int ngroups, int my_mpi_rank) {  /*group level; for now independently computed at current level*/
	for (int j = 0; j < ngroups; j++) {
		/* replacing each gcj−1 when j mod 32 = 0 using correct sck−1,k = j/32 */
		if ((j % block_size == 0) && (j != 0)) { /*the last one in the section*/
			gcj[j - 1] = sck[j / block_size - 1]; /* integer division, i.e. 10/3 = 3*/
		}

		if ((my_mpi_rank == 0) && (j == 0)) {
			gcj[j] = ggj[j] | (gpj[j] & c_minus_1);
		}
		else {
			gcj[j] = ggj[j] | (gpj[j] & gcj[j - 1]);
		}
	}


	return 0;
}

int c(int * gcj, int * ci, int * gi, int * pi, int bits, int my_mpi_rank){ /*lowest level; for now independently computed at current level*/

	for (int i = 0; i < bits; i++) {
		/* ci−1 = gcj−1 when i mod 32 = 0 and j = i/32 */
		if ((i % block_size == 0) && (i != 0)) { /*the last one in the section*/
			ci[i - 1] = gcj[i / block_size - 1]; /* integer division, i.e. 10/3 = 3*/
		}

		if ((my_mpi_rank == 0) && (i == 0)) {
			ci[i] = gi[i] | (pi[i] & c_minus_1);
		}
		else {
			ci[i] = gi[i] | (pi[i] & ci[i - 1]);
		}
	}

	return 0;

}

/*
 sum calculation
*/

int sum_cla(int * sumi, int * bin1, int * bin2, int * ci, int bits) { /*getting all the carries; do sum - (a XOR b) XOR c */
	for (int i = 0; i < bits; i++) {
		if (i == 0) {
			sumi[i] = (bin1[i] ^ bin2[i]) ^ c_minus_1;
		}
		else {
			sumi[i] = (bin1[i] ^ bin2[i]) ^ ci[i - 1];
		}
	}

	return 0;

}



