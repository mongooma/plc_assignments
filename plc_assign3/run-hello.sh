#!/bin/sh

srun --ntasks 256 --overcommit -o hello11111.log /gpfs/u/home/PCP8/PCP8mmnq/barn/a.out &
wait
