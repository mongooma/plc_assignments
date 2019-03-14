#!/bin/sh
#salloc --nodes 256 --time 60 --partition large # not necessary when not testing
#mpixlc -std=c99 -O main.c functions.c -o ./c.xl

sbatch --partition small --nodes 1 --time 5 ./sbatch.sh &
sbatch --partition small --nodes 2 --time 5 ./sbatch2.sh &
sbatch --partition small --nodes 4 --time 5 ./sbatch4.sh &
sbatch --partition small --nodes 8 --time 5 ./sbatch8.sh &
sbatch --partition small --nodes 16 --time 5 ./sbatch16.sh & 
sbatch --partition small --nodes 32 --time 5 ./sbatch32.sh &
sbatch --partition small --nodes 64 --time 5 ./sbatch64.sh &
sbatch --partition medium --nodes 128 --time 5 ./sbatch128.sh &
wait 
