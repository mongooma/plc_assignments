#!/bin/sh
srun --nodes=1 --ntasks=64  --overcommit -o result1.log ./c.xl &
wait
