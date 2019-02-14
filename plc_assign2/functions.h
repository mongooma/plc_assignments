#pragma once // check https://en.wikipedia.org/wiki/Include_guard for more info

void readInData(const char * filename,  
						char * hex_input_a, char * hex_input_b);

int convert_hex_2_bit(char * hex_input_a, char * hex_input_b, int * bin1, int bin2, int input_size);

int revert_binary(int * bin1, int * bin2, int bits);

int revert_hex_sum(int * sumi, int bits);

int convert_bit_2_hex(int * sumi, int bits);

int cla(bool use_barrier, int my_mpi_rank, int my_mpi_size, int alloc, int * bin1, int *bin2, int *sumi);

int g(int bits);

int p(int bits);

int gg(int ngroups);

int gp(int ngroups);

int sg(int nsections);

int sp(int nsections);

int ssg(int nsupersections);

int ssp(int supersection);

int ssc(int my_mpi_size, int my_mpi_rank, int nsupersections);

int sc(int nsections, int my_mpi_rank);

int gc(int ngroups, int my_mpi_rank);

int c(int bits, int my_mpi_rank);

int sum_cla(int bits);


