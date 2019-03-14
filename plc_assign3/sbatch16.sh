#!/bin/sh
srun --nodes=16 --ntasks=1024 --overcommit -o result16.log ./c.xl &
wait
