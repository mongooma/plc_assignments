#!/bin/sh
srun --nodes=128 --ntasks=8192 --overcommit -o result128.log ./c.xl &
wait
