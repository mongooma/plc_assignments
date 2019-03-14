#!/bin/sh


srun --nodes=5 --ntasks=320  --overcommit -o ./result1.log ./main.out --job-name=reduce_test -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com
