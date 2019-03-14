#!/bin/sh
srun --nodes=2 --ntasks=128 --overcommit -o result2.log ./c.xl &
wait
