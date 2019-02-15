#!/bin/bash

ranks=(2 4 8 16 32)
tests=(1 2 3 4 5 6 7) 
trials=(1 2 3 4 5 6 7 8 9 10)
use_barrier=(0 1)

printf "timestamp: `date`\n" >> log_cloud.txt

for trial in ${trials[*]};do 
	for rank in ${ranks[*]};do
		for t in ${tests[*]};do
			for b in ${use_barrier[*]};do
				printf "rank $rank, test $t, use_barrier $b" >> log_cloud.txt
				mpirun -np $rank ./a.out ./testdata/test\_input\_$t.txt $b >> log_cloud.txt
			done
		done
	done
done 
