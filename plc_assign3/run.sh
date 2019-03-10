mpicc -g -D DEBUG -Wall ./main.c ./functions.c
env RANK_NO=4 mpirun -np 4 ./a.out