#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=getAverageWig
#SBATCH --output=/scratch/las821/reports/slurm_getAverageWig_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_getAverageWig_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --mem=3GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu


val=$SLURM_ARRAY_TASK_ID
file=`sed -n ${val}p forAvg.txt`

echo "running command: perl /scratch/cgsb/ercan/scripts/chip/getAverageWig.pl $file"

perl /scratch/cgsb/ercan/scripts/chip/getAverageWig.pl $file

echo "done"

exit 0;

