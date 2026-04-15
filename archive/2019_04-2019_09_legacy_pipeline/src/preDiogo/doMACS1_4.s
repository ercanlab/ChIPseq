#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=MACS
#SBATCH --output=/scratch/las821/reports/slurm_MACS4_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_MACS4_%j.err
#SBATCH --time=05:00:00
#SBATCH --nodes=2
#SBATCH --mem=15GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load macs/1.4.3

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p forMacs.txt`
all=(`echo $params | sed 's/\s/ /g'`)

perl /scratch/cgsb/ercan/scripts/chip/slurm/macs_module_multiple_slurm.pl ${all[0]} ${all[1]} ${all[2]}
exit 0;
