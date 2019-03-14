#!/bin/sh
srun --nodes=8 --ntasks=512 --overcommit -o result8.log ./c.xl &
wait
