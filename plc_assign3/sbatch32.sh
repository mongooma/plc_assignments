#!/bin/sh
srun --nodes=32 --ntasks=2048 --overcommit -o result32.log ./c.xl &
wait
