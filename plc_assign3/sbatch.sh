#!/bin/sh

#!/bin/sh
#SBATCH --job-name=TESTING
#SBATCH -t 04:00:00
#SBATCH -D /gpfs/u/<home or barn or scratch>/<project>/<user>
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>

srun --nodes=1 --ntasks=64  --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
srun --nodes=2 --ntasks=128 --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
srun --nodes=4 --ntasks=256 --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
srun --nodes=8 --ntasks=512 --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
srun --nodes=16 --ntasks=1024 --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
srun --nodes=32 --ntasks=2048 --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
srun --nodes=64 --ntasks=4096 --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
srun --nodes=128 --ntasks=8192 --overcommit -o result.log ./main.xl --job-name=reduce_test -begin=01:00:00 -D --mail-type=ALL --mail-user=rhythm_qing@hotmail.com &
wait
# srun --ntasks  --overcommit -o hello.log /gpfs/u/home/SPNR/SPNRcaro/barn/mpi-hello.xl 