#!/bin/bash

ranks=(2 4 8 16 32)
tests=(1 2 3 4 5 6 7) 
trials=(1 2 3 4 5 6 7 8 9 10)

printf "timestamp: `date`\n" >> log_cloud.txt

for trial in ${trials[*]};do 
	for rank in ${ranks[*]};do
		for t in ${tests[*]};do
			printf "rank $rank, test $t, " >> log_cloud.txt
			mpirun -np $rank ./a.out ./testdata/test\_input\_$t.txt >> log_cloud.txt
		done
	done
done 
