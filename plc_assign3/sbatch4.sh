#!/bin/sh
srun --nodes=4 --ntasks=256 --overcommit -o result4.log ./c.xl &
wait
