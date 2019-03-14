#!/bin/sh

#mpixlc -O5 hello.c -o mpi-hello.xl

sbatch --partition debug --nodes 4 --time 5 ./run-hello.sh 

