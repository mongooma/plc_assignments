#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>
#include "functions.h"

// Add 1 to array size because strings must be null terminated
char hex_input_a[HEX_INPUT_SIZE+1]={0};
char hex_input_b[HEX_INPUT_SIZE+1]={0};

// EXAMPLE DATA STRUCTURE DESIGN AND LAYOUT FOR CLA
#define HEX_INPUT_SIZE 262144 /*  represent the 1,048,576 bit input numbers as 262,144 hex 
								digits*/ 
//#define input_size 1024
#define block_size 32

//Do not touch these defines 
#define digits (input_size+1) 
#define bits digits*4  /* every hex character is 4 digit binary */
#define ngroups bits/block_size 
#define nsections ngroups/block_size 
#define nsupersections nsections/block_size

//Global definitions of the various arrays used in steps for easy access 
int gi[bits] = {0}; 
int pi[bits] = {0}; 
int ci[bits] = {0}; /* what's this kind of initialization? */

int ggj[ngroups] = { 0 };
int gpj[ngroups] = { 0 }; 
int gcj[ngroups] = { 0 };

int sgk[nsections] = { 0 }; 
int spk[nsections] = { 0 }; 
int sck[nsections] = { 0 };

int ssgl[nsupersections] = { 0 }; 
int sspl[nsupersections] = { 0 }; 
int sscl[nsupersections] = { 0 };
int sumi[bits] = { 0 };

//Integer array of inputs in binary form 
int* bin1 = NULL; 
int* bin2 = NULL;

//Character array of inputs in hex form 
char* hex1 = NULL; 
char* hex2 = NULL;


// added declarations
int c_minus_1 = 0; /*primary value to set off the calculation chain*/
int i, j, k, l;
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

/*1st level - bits 4096 */
int g(int s, int e){  /*generates; 1 1*/
	
	for (i = s; i < e+1; i++) { 

		gi[i] = bin1[i] & bin2[i];
	}
	return 0;
}


int p(int s, int e) { /*propagates: 1 0, 0 1, 1 1*/
	for (i = s; i < e+1; i++) {

		pi[i] = bin1[i] | bin2[i];
	}
	return 0;
}

/*2nd level - group 512*/
int gg(int s, int e) {

	for (j = s; j < e+1; j++) {
		for(int b = 0+1; b < block_size; b++){
			int tmp = 1;
			for(int b_minus = 0; b_minus < b; b_minus++){ 
				if (b_minus == b-1) {tmp &= gi[j * block_size + (block_size - 1 - b_minus)]; /* b_minus maximum is block_size - 1*/}
				else{ tmp &= pi[j * block_size + (block_size - 1 - b_minus)];}
			}
			ggj[j] |= tmp;
		}
	}
	return 0;

}

int gp(int s, int e) { /*the logic will be the same for other propagate functions at higher levels:
			in order to propagate from previous group, the chain: 1-->2-->3-->4--> must be valid at every point
		   */
	for (j = s; j < e+1; j++) {
		gpj[j] = 1;
		for (int b = block_size-1; b > -1; b--){
			gpj[j] &= pi[j * block_size + b];
		}
	}
	return 0;

}

/* 3rd level - sections 64*/

int sg(int s, int e) { /*check gg()*/

	for (k = s; k < e+1; k++) {
		for(int b = 0+1; b < block_size; b++){
			int tmp = 1;
			for(int b_minus = 0; b_minus < b; b_minus++){ 
				if (b_minus == b-1) {tmp &= ggj[k * block_size + (block_size - 1 - b_minus)]; /* b_minus maximum is block_size - 1*/}
				else{ tmp &= gpj[k * block_size + (block_size - 1 - b_minus)];}
			}
			sgk[k] |= tmp;
		}
	}
	return 0;
}

int sp(int s, int e) {

	for (k = s; k < e+1; k++) {
		spk[k] = 1;
		for (int b = block_size-1; b > -1; b--){
			spk[k] &= gpj[k * block_size + b];
		}
	}
	return 0;

}

/* 4th level - supersection 8*/

int ssg(int s, int e) { /*check gg()*/
	for (l = s; l < e+1; l++) {
		for(int b = 0+1; b < block_size; b++){
			int tmp = 1;
			for(int b_minus = 0; b_minus < b; b_minus++){ 
				if (b_minus == b-1) {tmp &= sgk[l * block_size + (block_size - 1 - b_minus)]; /* b_minus maximum is block_size - 1*/}
				else{ tmp &= spk[l * block_size + (block_size - 1 - b_minus)];}
			}
			ssgl[l] |= tmp;
		}
	}
	return 0;
}

int ssp(int s, int e) { 

	for (l = s; l < e+1; l++) {
		sspl[l] = 1;
		for (int b = block_size-1; b > -1; b--){
			sspl[l] &= spk[l * block_size + b];
		}
	}
	return 0;

}

/* reduce: carries computation
*/

int ssc(int s, int e) { /*at supersection level*/ 
/* here we communicate between ranks since we need sscl[l-1] from neighbor node*/
/* Use MPI Isend and MPI Irecv routines to exchange the nearest-neighbor 
carry bit information at the sscl level*/
/* so that all processes share the information of ssc_{l-1}, l= k/block_size  
	for next step sc_{k-1} correction, k mod block_size = 0*/
/* send and recv could be called globally since source is stated*/
	
	for (l = s; l < e+1; l++) {
		if ((l == 0)) {
			sscl[l] = ssgl[l] | (sspl[l] & c_minus_1);
		}
		else {
			/* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
               int source, int tag, MPI_Comm comm, MPI_Request *request)*/
			MPI_Irecv(&(sscl[l-1]), 1，MPI_INT, l-1, /*message tag*/, MPI_COMM_WORLD, &status);
			sscl[l] = ssgl[l] | (sspl[l] & sscl[l-1]);

			/* int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
            int tag, MPI_Comm comm, MPI_Request *request)*/
			MPI_Isend(&(sscl[l]), 1, MPI_INT, l+1, /*message tag*/, ); 

			/* Question: do we need to broadcast the sscl variable?*/
			/* if not, then each rank will only maintain 2 blocks of sscl ? */
		}
	}

	return 0;

}

int sc() { /*section level; for now independently computed at current level*/

	for (k = 0; k < nsections; k++) {
		if (k == 0) {
			sck[k] = sgk[k] | (spk[k] & c_minus_1);
		}
		else {
			sck[k] = sgk[k] | (spk[k] & sck[k - 1]);
		}
		/*why? correct ssc_l, l = k / 8 for all sections k;*/
		if (k % block_size == block_size - 1) { /*the last one in the section*/
			sscl[k / block_size] = sck[k]; /* integer division, i.e. 10/3 = 3*/
		}
	}

	return 0;
}

int gc() {  /*group level; for now independently computed at current level*/
	for (j = 0; j < ngroups; j++) {
		if (j == 0) {
			gcj[j] = ggj[j] | (gpj[j] & c_minus_1);
		}
		else {
			gcj[j] = ggj[j] | (gpj[j] & gcj[j - 1]);
		}
		/*todo: why need to do this? correct sck_k, k = j / 8 for all sections j;*/
		if (j % block_size == block_size -1) { /*the last one in the section*/
			sck[j / block_size] = gcj[j]; /* integer division, i.e. 10/3 = 3*/
		}
	}


	return 0;
}

int c(){ /*lowest level; for now independently computed at current level*/
	for (i = 0; i < bits; i++) {
		if (i == 0) {
			ci[i] = gi[i] | (pi[i] & c_minus_1);
		}
		else {
			ci[i] = gi[i] | (pi[i] & ci[i - 1]);
		}
		/*correct gcj_j, j = i / 8 for all sections i; added another "or" condition for previous level*/
		if (i % block_size == block_size -1) { /*the last one in the section*/
			gcj[i / block_size ] = ci[i]; /* integer division, i.e. 10/3 = 3*/
		}
	}

	return 0;

}

/*
 sum calculation
*/

int sum_cla() { /*getting all the carries; do sum - (a XOR b) XOR c */
	for (i = 0; i < bits; i++) {
		if (i == 0) {
			sumi[i] = (bin1[i] ^ bin2[i]) ^ c_minus_1;
		}
		else {
			sumi[i] = (bin1[i] ^ bin2[i]) ^ ci[i - 1];
		}
	}

	return 0;

}

/* master */
/* cla executed within one rank */
int cla(bool use_barrier, int my_mpi_rank, int my_mpi_size, int alloc){

	/* all ranks except Rank 0, should post an MPI Irecv message before the cla 
	calculations start */ 
	/* what inputs are they waiting for ?*/
	/* is this the c-minus for the left most group/section/nsection ?*/

	MPI_Request status;
	int c_minus_1 = 0;

	if(my_mpi_rank != 0){

		/*  int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
               int source, int tag, MPI_Comm comm, MPI_Request *request)*/
		MPI_Irecv(/*buf ?*/, 1, /*dtype*/, /*source*/, /*tag*/, MPI_COMM_WORLD, &status);
		/* call MPI_Wait, etc to check the status */
	}	
	

	/* perform cla for each rank*/

	/* g */
	/* p */
	g(my_mpi_rank * alloc, (my_mpi_rank + 1) * alloc - 1);
	p(my_mpi_rank * alloc, (my_mpi_rank + 1) * alloc - 1);

	/* between each algorithm step (level) perform an MPI Barrier collective operation to 
	keep all ranks in step with each other. */

	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

	/* gg */
	/* gp */
	gg(my_mpi_rank * (ngroups / my_mpi_size) /* start group*/, (my_mpi_rank + 1) *(ngroups / my_mpi_size) - 1 /* end group*/);
	gp(my_mpi_rank * (ngroups / my_mpi_size), (my_mpi_rank + 1) *(ngroups / my_mpi_size) - 1)

	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

	/* sg */
	/* sp */
	sg(my_mpi_rank * (nsections / my_mpi_size), (my_mpi_rank + 1) *(nsections / my_mpi_size) - 1);
	sp(my_mpi_rank * (nsections / my_mpi_size), (my_mpi_rank + 1) *(nsections / my_mpi_size) - 1);

	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

	/* ssg */
	/* ssp */
	ssg(my_mpi_rank * (nsupersections / my_mpi_size), (my_mpi_rank + 1) *(nsupersections / my_mpi_size) - 1);
	ssp(my_mpi_rank * (nsupersections / my_mpi_size), (my_mpi_rank + 1) *(nsupersections / my_mpi_size) - 1);


	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

	/* ssc */
	ssc(my_mpi_rank * (nsupersections / my_mpi_size), (my_mpi_rank + 1) *(nsupersections / my_mpi_size) - 1);


	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

	/* sc */
	sc(my_mpi_rank * (nsections / my_mpi_size), (my_mpi_rank + 1) *(nsections / my_mpi_size) - 1);


	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

	/* gc */

	gc(my_mpi_rank * (ngroups / my_mpi_size), (my_mpi_rank + 1) *(ngroups / my_mpi_size) - 1)

	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

	/* c */

	c(my_mpi_rank * alloc, (my_mpi_rank + 1) * alloc - 1);

	if(use_barrier){MPI_Barrier(); /* make it optional for performance study*/}

}


int convert_bit_2_hex(){
	/* implement a look-up table instead */
	
	char *tmp1, *tmp2, tmp_[2];
	tmp_[1] = '\0';

	bin1 = calloc(bits, sizeof(int));
	bin2 = calloc(bits, sizeof(int));

	/* hex1 and hex2 of length input_size + 1*/
	for (i = 0; i < input_size + 1; i++) {
		tmp_[0] = hex1[i];
		tmp1 = lookup[strtol(tmp_, 0, 16)];
		tmp_[0] = hex2[i];
		tmp2 = lookup[strtol(tmp_, 0, 16)];
		for (j = 0; j < 4; j++) {
			tmp_[0] = tmp1[j];
			bin1[i * 4 + j] = atoi(tmp_); /* atoi will return exact as it looks*/
			tmp_[0] = tmp2[j];
			bin2[i * 4 + j] = atoi(tmp_); /* atoi will return exact as it looks*/
		}
	}

}


int revert_binary(){
	/* revert the binary string to set the last digit as the MSB*/
	for (i = 0; i < (bits / 2); i++) {

		int tmp = bin1[i];
		bin1[i] = bin1[(bits - 1) - i];
		bin1[(bits - 1) - i] = tmp;
		tmp = bin2[i];
		bin2[i] = bin2[(bits - 1) - i];
		bin2[(bits - 1) - i] = tmp;
	}

}

int revert_hex_sum(){

	for (i = 0; i < bits / 2; i++) {
		int tmp = sumi[i];
		sumi[i] = sumi[(bits - 1) - i];
		sumi[(bits - 1) - i] = tmp;
	}

}

int convert_hex_2_bit(){

	for (i = 0; i < bits; i += 4) {
		int tmp = 0; 
		for (j = 0; j < 4; j++) {
			tmp += (sumi[i + j] % 2);
			if (j == 3) break;
			tmp = tmp * 2;
		}
	}


}