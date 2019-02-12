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

// EXAMPLE DATA STRUCTURE DESIGN AND LAYOUT FOR CLA 
#define input_size 8

//Do not touch these defines 
#define digits (input_size+1) 
#define bits digits*4  /* every hex character is 4 digit binary */

//Global definitions of the various arrays used in steps for easy access 
int gi[bits] = { 0 };
int pi[bits] = { 0 };
int ci[bits] = { 0 }; /* what's this kind of initialization? */

int sumi[bits] = { 0 };

//Integer array of inputs in binary form 
int* bin1 = NULL;
int* bin2 = NULL;

//Character array of inputs in hex form 
char* hex1 = NULL;
char* hex2 = NULL;


// added declarations
int c_minus_1 = 0; /*primary value to set off the calculation chain*/
int i, j;
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
int g() { /*todo: check input form; might need to perform this 4 in a group*/

	for (i = 0; i < bits; i++) {
		/*change hex to binary*/

		gi[i] = bin1[i] & bin2[i];
	}
	return 0;
}
int p() {
	for (i = 0; i < bits; i++) {
		/*change hex to binary*/

		pi[i] = bin1[i] | bin2[i];
	}
	return 0;
}


int c() {
	for (i = 0; i < bits; i++) {
		if (i == 0) {
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

int sum_cla() {
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

int cla() { /* todo: check return value*/

	g();
	p();
	c();
	sum_cla();

	return 0;

}

int main() {

	//1. Read two 1024 digit hex numbers as input. Add a leading hex “0” digit to input. /*todo: why?*/
	// ./cla < test0.txt


	hex1 = calloc(digits + 2, sizeof(char)); /* + 1 for '\0' */
	hex2 = calloc(digits + 2, sizeof(char));
	hex1[0] = '0'; /*pend a leading 0, as MSB (the most significant bit, highest value position)*/
	hex2[0] = '0';

	scanf("%s\n%s\n", &(*(hex1 + 1)), &(*(hex2 + 1))); /* could also use &(hex1[1]) here */ /*input from stdin*/
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

	cla();

	//6. reconvert the binary sum of the two numbers into hex output.

#ifdef DEBUG

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

#ifdef DEBUG_pass

	for (i = 0; i < bits; i++) {
		printf("%d", ci[i]);
	}
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
		int tmp=0;
		for (j = 0; j < 4; j++) {
			tmp += (sumi[i + j] % 2);
			if (j == 3) break;
			tmp = tmp * 2;
		}
		printf("%x", tmp); /* %x format for hex*/
	}
	printf("\n");


	free(hex1);
	free(hex2);
	free(bin1);
	free(bin2);


	return 0;
}