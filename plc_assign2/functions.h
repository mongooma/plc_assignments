#pragma once // check https://en.wikipedia.org/wiki/Include_guard for more info

void readInData(const char * filename,  
						char * hex_input_a, char * hex_input_b);

int convert_hex_2_bit(char * hex_input_a, char * hex_input_b, int * bin1, int *bin2, int input_size);

int revert_binary(int * bin1, int * bin2, int bits);

int revert_hex_sum(int * sumi, int bits);

int convert_bit_2_hex(int * sumi, int bits);

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


