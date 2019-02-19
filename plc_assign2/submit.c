#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<mpi.h>

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


#define HEX_INPUT_SIZE 262144
// #define HEX_INPUT_SIZE 64

/* variables for performance report:
	ranks: 2, 4, 8, 16, 32
	use_barrier: 0, 1

*/

// Compile Code: mpicc -g -Wall main.c 
// Example Run Code: mpirun -np 4 ./a.out input.txt


// Both input and output files will 524490 bytes in size. The two additional
//    characters are due to a newline characters between each input and at the
//    end of the file.



/*
 use 32 bit blocks and extend it to work in parallel for a 1M (e.g., 1,048,576) 
 bit CLA adder with 32 bit blocks using up to 32 MPI ranks on the class server, 
 mastiff.cs.rpi.edu. 
*/


/*notation: significance order: b<--|<--a*/


void readInData(const char * filename,  
						char * hex_input_a, char * hex_input_b);

int convert_hex_2_bit(char * hex_input_a, char * hex_input_b, int * bin1, int *bin2, int input_size);

int revert_binary(int * bin1, int * bin2, int bits);

int revert_hex_sum(int * sumi, int bits);

int convert_bit_2_hex(int * sumi, int bits, char * output);

int cla(int use_barrier, int my_mpi_rank, int my_mpi_size, int alloc, int * bin1, int *bin2, int *sumi);

int g(int bits, int *gi, int *bin1, int * bin2);

int p(int bits, int *pi, int *bin1, int * bin2);

int gg(int ngroups, int *ggj, int *gi, int * pi);

int gp(int ngroups, int *gpj, int * pi);

int sg(int nsections, int *sgk, int *ggj, int * gpj);

int sp(int nsections, int *spk, int *gpj);

int ssg(int nsupersections, int *ssgl, int *sgk, int * spk);

int ssp(int nsupersections, int *sspl, int * spk);

int ssc(int * sscl, int * ssgl, int * sspl, int my_mpi_size, int my_mpi_rank, int nsupersections);

int sc(int * sscl, int * sck, int * sgk, int * spk, int nsections, int my_mpi_rank);

int gc(int * sck, int * gcj, int * ggj, int * gpj, int ngroups, int my_mpi_rank);

int c(int * gcj, int * ci, int * gi, int * pi, int bits, int my_mpi_rank);

int sum_cla(int * sumi, int * bin1, int * bin2, int * ci, int bits);


/*notation: significance order: b<--|<--a*/




void readInData(const char * filename,  
						char * hex_input_a, char * hex_input_b){
	/* to debug */

	FILE *my_input_file=NULL;

	#ifdef DEBUG_MUTE_PASS
	#endif

	if( (my_input_file = fopen( filename, "r")) == NULL ){
		printf("Failed to open input data file: %s \n", filename);
		exit(-1);
	}

	fscanf( my_input_file, "%s %s", hex_input_a, hex_input_b );

	#ifdef DEBUG_MUTE_PASS
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

int convert_bit_2_hex(int * sumi, int bits, char * output){

	FILE *f = fopen(output, "w"); /* use fprintf for debugging output*/
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
	int * ci = calloc(bits + 1, sizeof(int));
	ci += 1; //this allows for ci[-1] indexing  

	int * ggj = calloc(ngroups, sizeof(int));
	int * gpj = calloc(ngroups, sizeof(int));
	int * gcj = calloc(ngroups + 1, sizeof(int));
	gcj += 1;

	int * sgk = calloc(nsections, sizeof(int));
	int * spk = calloc(nsections, sizeof(int));
	int * sck = calloc(nsections + 1, sizeof(int));
	sck += 1;

	int * ssgl = calloc(nsupersections, sizeof(int));
	int * sspl = calloc(nsupersections, sizeof(int));
	int * sscl = calloc(nsupersections + 1, sizeof(int));
	sscl += 1;

	#ifdef DEBUG_MUTE_PASS_p
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
		#ifdef DEBUG 
		if (i < bits){
			printf("g: gi[%d]=%d \n", i, gi[i]);
		}
		#endif
	}

	return 0;
}


int p(int bits, int *pi, int *bin1, int * bin2) { /*propagates: 1 0, 0 1, 1 1*/
	for (int i = 0; i < bits; i++) {

		pi[i] = bin1[i] | bin2[i];
		#ifdef DEBUG 
		if (i < bits){
			printf("p: pi[%d]=%d \n", i, pi[i]);
		}
		#endif

	}
	return 0;
}

/*2nd level - group 512*/
int gg(int ngroups, int *ggj, int *gi, int * pi) {

	for (int j = 0; j < ngroups; j++) {
		for(int b = 0+1; b < block_size+1; b++){ //blocksize+1
			int tmp = 1;
			for(int b_minus = 0; b_minus < b; b_minus++){ 
				if (b_minus == b-1) { tmp &= gi[j * block_size + (block_size - 1 - b_minus)];
											/* b_minus maximum is block_size - 1*/}
								else{ tmp &= pi[j * block_size + (block_size - 1 - b_minus)];
								}
			}
			ggj[j] |= tmp;
		}
		#ifdef DEBUG_MUTE_PASS
		if(j == 0){
			printf("gg: gg[0] = %d\n", ggj[0]);
		}
		#endif
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
		#ifdef DEBUG_MUTE_PASS
		if(j == 0){
			printf("gp: gp[0] = %d\n", gpj[0]);
		}
		#endif
	}
	return 0;

}

/* 3rd level - sections 64*/

int sg(int nsections, int *sgk, int *ggj, int * gpj) { /*check gg()*/

	for (int k = 0; k < nsections; k++) {
		for(int b = 0+1; b < block_size+1; b++){
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
		for(int b = 0+1; b < block_size+1; b++){
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

	#ifdef DEBUG_MUTE_PASS_p
		printf("rank %d, ssc: here0.\n", my_mpi_rank);
	#endif


	for(int i = 0; i < my_mpi_size; i++){

		if ((my_mpi_rank == i) && ((i == 0) || (my_mpi_size == 1)) ){ 

			sscl[-1] = c_minus_1;

			for (int l = 0; l < nsupersections; l++) {

				sscl[l] = ssgl[l] | (sspl[l] & sscl[l-1]);
				#ifdef DEBUG_MUTE
					printf("ssc: actual sscl, mpi=1, sscl[%d] = %d \n", l, sscl[l]);
				#endif

			}

			if(my_mpi_size != 1){
				MPI_Isend(sscl + nsupersections - 1, 1, MPI_INT, my_mpi_rank + 1, 0, MPI_COMM_WORLD, &request); 
			}else{ return 0;}

		} 
		
		MPI_Barrier(MPI_COMM_WORLD); /*To have every process synchronized*/ 


		if( (my_mpi_rank == i) && ( i > 0) ){
			/* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
           int source, int tag, MPI_Comm comm, MPI_Request *request)*/
			MPI_Irecv(sscl - 1, 1, MPI_INT, my_mpi_rank - 1, 0, MPI_COMM_WORLD, &request);
			// update buf as initial carrier for this rank (rank > 0)

			#ifdef DEBUG
			printf("rank %d, ssc: recv sscl[-1] = %d \n", my_mpi_rank, sscl[-1]);
			#endif

			for (int l = 0; l < nsupersections; l++) {
				sscl[l] = ssgl[l] | (sspl[l] & sscl[l-1]); 
			}

		}

		

		MPI_Barrier(MPI_COMM_WORLD); /*To have every process synchronized*/ 

		
		if((my_mpi_rank == i) && (i < my_mpi_size-1)){
		
			/* int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
			int tag, MPI_Comm comm, MPI_Request *request)*/
			MPI_Isend(sscl + nsupersections - 1, 1, MPI_INT, my_mpi_rank+1, 0, MPI_COMM_WORLD, &request); 
			#ifdef DEBUG
				printf("rank %d, ssc: send sscl[end] = %d \n", my_mpi_rank, sscl[nsupersections-1]);
			#endif
		}

		MPI_Barrier(MPI_COMM_WORLD); /*To have every process synchronized*/ 


	}

	return 0;

}

int sc(int *sscl, int * sck, int * sgk, int * spk, int nsections, int my_mpi_rank) { /*section level; for now independently computed at current level*/

	for (int k = 0; k < nsections; k++) {

		// /* replacing each sck−1 when k mod 32 = 0 using correct s2cl−1,l = k div 32 */
		if (k % block_size == 0) { /* get carrier from last sscl*/
			
			sck[k-1] = sscl[k / block_size - 1];

			#ifdef DEBUG
			printf("rank: %d, sc: sck[-1] = %d \n", my_mpi_rank, sck[k-1]);
			#endif

		}

		sck[k] = sgk[k] | (spk[k] & sck[k-1]);//


	}

	return 0;
}


int gc(int * sck, int * gcj, int * ggj, int * gpj, int ngroups, int my_mpi_rank) {  /*group level; for now independently computed at current level*/
	for (int j = 0; j < ngroups; j++) {
		// /* replacing each gcj−1 when j mod 32 = 0 using correct sck−1,k = j/32 */
		if (j % block_size == 0) { /*the last one in the section*/
			
			gcj[j - 1] = sck[j / block_size - 1]; /* integer division, i.e. 10/3 = 3*/
			
			#ifdef DEBUG
			printf("rank: %d, gc: gcj[-1] = %d \n", my_mpi_rank, gcj[j-1]);
			#endif
		}
		
		gcj[j] = ggj[j] | (gpj[j] & gcj[j - 1]);
		
	}


	return 0;
}

int c(int * gcj, int * ci, int * gi, int * pi, int bits, int my_mpi_rank){ /*lowest level; for now independently computed at current level*/

	for (int i = 0; i < bits; i++) {
		/* ci−1 = gcj−1 when i mod 32 = 0 and j = i/32 */
		if (i % block_size == 0){ /*the last one in the section*/
			
			ci[i - 1] = gcj[i / block_size - 1]; /* integer division, i.e. 10/3 = 3*/

			#ifdef DEBUG
			printf("rank: %d, c: ci[-1] = %d \n", my_mpi_rank, ci[i-1]);
			#endif 
		}

		ci[i] = gi[i] | (pi[i] & ci[i - 1]);

	}

	return 0;

}

/*
 sum calculation
*/

int sum_cla(int * sumi, int * bin1, int * bin2, int * ci, int bits) { /*getting all the carries; do sum - (a XOR b) XOR c */

	for (int i = 0; i < bits; i++) {
		
		sumi[i] = (bin1[i] ^ bin2[i]) ^ ci[i - 1];
		
	}

	return 0;

}





int main(int argc, char ** argv){
	/*
	if( argc != 2 ){
		printf("Not sufficient arguments, only %d found. \n", argc);
		exit(-1);
	}
	*/
	// Add 1 to array size because strings must be null terminated
	double start_time, end_time;

	int use_barrier = atoi(argv[argc-1]);
	//printf("Use use_barrier: %d", use_barrier);

	char * hex_input_a = calloc(HEX_INPUT_SIZE+1, sizeof(char)); //add 1 for '\0'
	char * hex_input_b = calloc(HEX_INPUT_SIZE+1, sizeof(char));

	//Integer array of inputs in binary form 4;
	int * bin1 = calloc((HEX_INPUT_SIZE) * 4, sizeof(int)); /* keep the fashion of assign 1, 0000 for the */
	int * bin2 = calloc((HEX_INPUT_SIZE) * 4, sizeof(int));


	int my_mpi_size = -1; // total no. of ranks
	int my_mpi_rank = -1; // the current rank of the process

	MPI_Init( &argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

	int allocation = (HEX_INPUT_SIZE * 4) / my_mpi_size;
	int * bin_rank_1 = calloc(allocation, sizeof(int));
	int * bin_rank_2 = calloc(allocation, sizeof(int));
	int * sumi = calloc(allocation, sizeof(int));
	int * sumi_all = calloc(HEX_INPUT_SIZE * 4 , sizeof(int));

	/*1. MPI rank 0 read in bit number, */

	if(my_mpi_rank == 0){

		readInData(argv[1], hex_input_a, hex_input_b); /* not debugged */

		/*2. MPI rank 0 convert hex to binary number, revert,  */

		convert_hex_2_bit(hex_input_a, hex_input_b, bin1, bin2, HEX_INPUT_SIZE); /* not debugged */

		revert_binary(bin1, bin2, HEX_INPUT_SIZE * 4);	

	}

	#ifdef DEBUG
	if(my_mpi_rank == 3){
		for(int i=0; i < 4; i ++){
			printf("bin1[%d]: %d \n", i, sumi[i]);
		}
	}
	#endif


	/*3. MPI Rank 0 distribute the input binary arrays to each rank in the correct order 
	where (for a 32 ranks conﬁguration) MPI rank 1 has bits 32,768 through 65,535 
	and MPI rank 2 has bits 65536 through 98304 and so on. */

	/* int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
        void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
        MPI_Comm comm)*/
	#ifdef DEBUG_p
		printf("rank %d, main: here00! \n", my_mpi_rank); /* three reached here */
	#endif 
	MPI_Scatter(bin1, allocation, MPI_INT, bin_rank_1, allocation, MPI_INT, 0, MPI_COMM_WORLD);
	#ifdef DEBUG_p
		printf("rank %d, main: here00! \n", my_mpi_rank); /* three reached here */
	#endif 
	MPI_Scatter(bin2, allocation, MPI_INT, bin_rank_2, allocation, MPI_INT, 0, MPI_COMM_WORLD);
	#ifdef DEBUG_p
		printf("rank %d, main: here01! \n", my_mpi_rank); /* three reached here */
	#endif 

	free(hex_input_a);
	free(hex_input_b);
	free(bin1);
	free(bin2);


	#ifdef DEBUG_P
		printf("rank %d, main: here1! \n", my_mpi_rank); /* three reached here */
	#endif 

	MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/



	/*execute algorithm -> see cla() */

	start_time = MPI_Wtime(); /* a time tick fashion*/

	cla(use_barrier, my_mpi_rank, my_mpi_size, allocation, bin_rank_1, bin_rank_2, sumi); 

	end_time = MPI_Wtime();

	#ifdef DEBUG_P
		printf("rank %d, main: here2! \n", my_mpi_rank); /* three reached here */
	#endif 


	/* for synchronization */
	MPI_Barrier(MPI_COMM_WORLD); /* make it optional for performance study*/


	/* Have each rank send their part of the ﬁnal sumi solution to Rank 0. */

	/*  MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm) */	
	
	#ifdef DEBUG
	if(my_mpi_rank == 3){
		for(int i=0; i < 4; i ++){
			printf("sumi[%d]: %d \n", i, sumi[i]);
		}
	}
	#endif

	MPI_Gather(sumi, allocation, MPI_INT, sumi_all, allocation, MPI_INT, 0, MPI_COMM_WORLD); 
	

	/*  Rank 0 will then re-reverse the sum and output the ﬁnal result. */
	if( my_mpi_rank == 0){
		revert_hex_sum(sumi_all, HEX_INPUT_SIZE * 4);
		convert_bit_2_hex(sumi_all, HEX_INPUT_SIZE * 4, argv[2]); // submit otuput file
		end_time = MPI_Wtime();
	}


	MPI_Finalize();


	// if(my_mpi_rank == 0){ /*only monitor root process rank 0*/
	// 	printf("time in seconds: %f, (use_barrier: %d, ranks: %s)\n", 
	// 							(end_time - start_time), use_barrier, argv[0] );
	// }

	free(bin_rank_1);
	free(bin_rank_2);
	free(sumi);
	free(sumi_all);


}