/* Carry-Lookahead Adder

simulates a 4,096 bit Carry Lookahead
Adder using 8 bit blocks. 


notice at this point no parallelization technique is used.
(higher levels, g, s, ss are not used for now, 
thus for all the same level computations, use one for loop and set initial carry as 0)

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// EXAMPLE DATA STRUCTURE DESIGN AND LAYOUT FOR CLA 
#define input_size 262144
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
int g(){  /*generates; 1 1*/
	
	for (i = 0; i < bits; i++) { 

		gi[i] = bin1[i] & bin2[i];
	}
	return 0;
}


int p() { /*propagates: 1 0, 0 1, 1 1*/
	for (i = 0; i < bits; i++) {

		pi[i] = bin1[i] | bin2[i];
	}
	return 0;
}

/*2nd level - group 512*/
int gg() {
	for (j = 0; j < ngroups; j++) {
		/* the logic will be the same for other g functions in higher-levels:
			generate (create a carry independently) -> propagate (create a carry if getting a carry from upstream)

			for a chain of 4 (1 -> 2 -> 3 -> 4), one possiblity to have a carry at the end of the chain is:
				1 -> 2(generate) -> 3(propagate) -> 4 (propagate) --> END
				
			*/
		ggj[j] = gi[j * block_size + 7] | 
				(pi[j * block_size + 7] & gi[j * block_size + 6]) |
				(pi[j * block_size + 7] & pi[j * block_size + 6] & gi[j * block_size + 5]) |
				(pi[j * block_size + 7] & pi[j * block_size + 6] & pi[j * block_size + 5] & pi[j * block_size + 4]) |
				(pi[j * block_size + 7] & pi[j * block_size + 6] & pi[j * block_size + 5] & pi[j * block_size + 4] & gi[j * block_size + 3]) |
				(pi[j * block_size + 7] & pi[j * block_size + 6] & pi[j * block_size + 5] & pi[j * block_size + 4] & pi[j * block_size + 3] & gi[j * block_size + 2]) |
				(pi[j * block_size + 7] & pi[j * block_size + 6] & pi[j * block_size + 5] & pi[j * block_size + 4] & pi[j * block_size + 3] & pi[j * block_size + 2] & gi[j * block_size + 1]) |
				(pi[j * block_size + 7] & pi[j * block_size + 6] & pi[j * block_size + 5] & pi[j * block_size + 4] & pi[j * block_size + 3] & pi[j * block_size + 2] & pi[j * block_size + 1] & gi[j * block_size]) 
			;
	}
	return 0;
}

int gp() { /*the logic will be the same for other propagate functions at higher levels:
			in order to propagate from previous group, the chain: 1-->2-->3-->4--> must be valid at every point
		   */
	for (j = 0; j < ngroups; j++) {
		gpj[j] = pi[j * block_size + 7] & pi[j * block_size + 6] & pi[j * block_size + 5] & pi[j * block_size + 4] & pi[j * block_size + 3] & pi[j * block_size + 2] & pi[j * block_size + 1] & pi[j * block_size + 0];
	}
	return 0;

}

/* 3rd level - sections 64*/

int sg() { /*check gg()*/
	for (k = 0; k < nsections; k++) {
		sgk[k] = ggj[k * block_size + 7] |
			(gpj[k * block_size + 7] & ggj[k * block_size + 6]) |
			(gpj[k * block_size + 7] & gpj[k * block_size + 6] & ggj[k * block_size + 5]) |
			(gpj[k * block_size + 7] & gpj[k * block_size + 6] & gpj[k * block_size + 5] & gpj[k * block_size + 4]) |
			(gpj[k * block_size + 7] & gpj[k * block_size + 6] & gpj[k * block_size + 5] & gpj[k * block_size + 4] & ggj[k * block_size + 3]) |
			(gpj[k * block_size + 7] & gpj[k * block_size + 6] & gpj[k * block_size + 5] & gpj[k * block_size + 4] & gpj[k * block_size + 3] & ggj[k * block_size + 2]) |
			(gpj[k * block_size + 7] & gpj[k * block_size + 6] & gpj[k * block_size + 5] & gpj[k * block_size + 4] & gpj[k * block_size + 3] & gpj[k * block_size + 2] & ggj[k * block_size + 1]) |
			(gpj[k * block_size + 7] & gpj[k * block_size + 6] & gpj[k * block_size + 5] & gpj[k * block_size + 4] & gpj[k * block_size + 3] & gpj[k * block_size + 2] & gpj[k * block_size + 1] & ggj[k * block_size])
			;
	}
	return 0;
}

int sp() { /*check gp()*/
	for (k = 0; k < nsections; k++) {
		spk[k] = gpj[k * block_size + 7] & gpj[k * block_size + 6] & gpj[k * block_size + 5] & gpj[k * block_size + 4] & gpj[k * block_size + 3] & gpj[k * block_size + 2] & gpj[k * block_size + 1] & gpj[k * block_size + 0];
	}
	return 0;

}

/* 4th level - supersection 8*/

int ssg() { /*check gg()*/
	for (l = 0; l < nsections; l++) {
		ssgl[l] = sgk[l * block_size + 7] |
			(spk[l * block_size + 7] & sgk[l * block_size + 6]) |
			(spk[l * block_size + 7] & spk[l * block_size + 6] & sgk[l * block_size + 5]) |
			(spk[l * block_size + 7] & spk[l * block_size + 6] & spk[l * block_size + 5] & spk[l * block_size + 4]) |
			(spk[l * block_size + 7] & spk[l * block_size + 6] & spk[l * block_size + 5] & spk[l * block_size + 4] & sgk[l * block_size + 3]) |
			(spk[l * block_size + 7] & spk[l * block_size + 6] & spk[l * block_size + 5] & spk[l * block_size + 4] & spk[l * block_size + 3] & sgk[l * block_size + 2]) |
			(spk[l * block_size + 7] & spk[l * block_size + 6] & spk[l * block_size + 5] & spk[l * block_size + 4] & spk[l * block_size + 3] & spk[l * block_size + 2] & sgk[l * block_size + 1]) |
			(spk[l * block_size + 7] & spk[l * block_size + 6] & spk[l * block_size + 5] & spk[l * block_size + 4] & spk[l * block_size + 3] & spk[l * block_size + 2] & spk[l * block_size + 1] & sgk[l * block_size])
			;
	}
	return 0;
}

int ssp() { /*check gp()*/
	for (l = 0; l < nsections; l++) {
		sspl[l] = spk[l * block_size + 7] & spk[l * block_size + 6] & spk[l * block_size + 5] & spk[l * block_size + 4] & spk[l * block_size + 3] & spk[l * block_size + 2] & spk[l * block_size + 1] & spk[l * block_size + 0];
	}
	return 0;

}

/* reduce: carries computation
*/

int ssc() { /*at supersection level*/
	
	for (l = 0; l < nsupersections; l++) {
		if (l == 0) {
			sscl[l] = ssgl[l] | (sspl[l] & c_minus_1);
		}
		else {
			sscl[l] = ssgl[l] | (sspl[l] & sscl[l-1]);
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
		if (k % 8 == 7) { /*the last one in the section*/
			sscl[k / 8] = sck[k]; /* integer division, i.e. 10/3 = 3*/
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
		if (j % 8 == 7) { /*the last one in the section*/
			sck[j / 8] = gcj[j]; /* integer division, i.e. 10/3 = 3*/
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
		if (i % 8 == 7) { /*the last one in the section*/
			gcj[i / 8] = ci[i]; /* integer division, i.e. 10/3 = 3*/
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

/* master*/

int cla() { /*master routine*/

	g();
	p();
	//gg();
	//gp();
	//sg();
	//sp();
	//sc();
	//gc();
	c();
	sum_cla();

	return 0;

}

int main(int argc, char ** argv) {

	clock_t start_t, end_t, total_t;

	setvbuf(stdout, NULL, _IONBF, 0);

	//1. Read two 1024 digit hex numbers as input. Add a leading hex “0” digit to input. /*todo: why?*/
	// ./cla < test0.txt


	hex1 = calloc(digits + 2, sizeof(char)); /* + 1 for '\0' */
	hex2 = calloc(digits + 2, sizeof(char));
	hex1[0] = '0'; /*pend a leading 0, as MSB (the most significant bit, highest value position)*/
	hex2[0] = '0';

	FILE *my_input_file=fopen(argv[1], "r");

	fscanf(my_input_file, "%s\n%s\n", &(*(hex1 + 1)), &(*(hex2 + 1))); /* could also use &(hex1[1]) here */ /*input from stdin*/
#ifdef DEBUG_pass
	printf("%s \n", hex2);
#endif // DEBUG

	/*2. Convert the “hex”(base - 16) string to a “binary”(base - 2) string.NOTE: That hex[0] is the most signiﬁcant digit, so when you convert to the “binary” string, you need to make sure that binary[4095] is the most signiﬁcant digit.This can be done by translating each hex digit into it’s binary format(e.g., B = 1011) and then reversing the whole binary array to make sure that binary[4095] is the most signiﬁcant digit and in the right order with respect to how the hex digits where translated.
		notice: different from normal hex-bin conversion, in which:

		int v = strtol(hex, 0, 16); // convert the hex value to a number

		for (i = 0; i < bits; i++) {
			bin1[i] = (v % 2 == 1) ? "1" : "0";
			v -= (v % 2);
			v /= 2;
			i++;
		}

		here the conversion should be:
		03A -> revert(0000 0011 1010) -> 1010 0011 0000

		implement a look-up table instead
	*/


	{ /* adding a mini scope here to destroy some tmp variables; */
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

#ifdef DEBUG_pass
	
	printf("%c \n", hex2[input_size]);
	for (i = 0; i < bits; i++) {
		printf("%d", bin2[i]);
}
	printf("\n");

	for (i = 0; i < bits; i++) {
		printf("%d", bin1[i]);
	}
	printf("\n");
#endif // DEBUG

	
	/* revert the binary string to set the last digit as the MSB*/
	for (i = 0; i < (bits / 2); i++) {

		int tmp = bin1[i];
		bin1[i] = bin1[(bits - 1) - i];
		bin1[(bits - 1) - i] = tmp;
		tmp = bin2[i];
		bin2[i] = bin2[(bits - 1) - i];
		bin2[(bits - 1) - i] = tmp;
	}


#ifdef DEBUG_pass
	printf("%d \n", bin2[0]);
#endif // DEBUG

	//3. Each bit can be a full integer and so you’ll have arrays of unsigned integers for each of the components parts of the adder.
	//4. functions for each step in the above algorithm. You can do it using for loops and do not have to perform the equation substitutions by hand as we did in class.
	//5. A master “cla” routine will run thru all the steps.
	
	start_t = clock();

	cla();

	end_t = clock();

	//6. reconvert the binary sum of the two numbers into hex output.

#ifdef DEBUG
	printf("\n");
	for (i = 0; i < bits; i++) {
		printf("%d", sumi[i]);
	}
	printf("\n");
#endif // DEBUG

	for (i = 0; i < bits / 2; i++) {
		int tmp = sumi[i];
		sumi[i] = sumi[(bits - 1) - i];
		sumi[(bits - 1) - i] = tmp;
	}
	
	//you can check your answer by creating a simple ripple carry adder based on computing the ci = gi + pi ∗ci−1 for all i and then computing the ﬁnal sum, based on sumi = ai ⊕bi ⊕ci−1. 


	//7. Convert the binary sumi array to “hex” and output your result in 1024 digit “hex” format.You don’t need to worry about overﬂow conditions. Note, you’ll need to re - invert / reverse the answer so that it prints correctly.

	for (i = 0; i < bits; i += 4) {
		int tmp = 0; 
		for (j = 0; j < 4; j++) {
			tmp += (sumi[i + j] % 2);
			if (j == 3) break;
			tmp = tmp * 2;
		}
		//printf("%X", tmp); /* %x format for hex*/
	}
	//printf("\n");

	printf("time in seconds: %f \n", 
								(double)(end_t - start_t)/CLOCKS_PER_SEC);

	free(hex1);
	free(hex2);
	free(bin1);
	free(bin2);


	return 0;
}