#!/bin/sh
srun --nodes=64 --ntasks=4096 --overcommit -o result64.log ./c.xl &
wait
