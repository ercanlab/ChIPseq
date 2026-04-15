#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=normWig
#SBATCH --output=/scratch/las821/reports/slurm_normWig_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_normWig_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=2
#SBATCH --mem=5GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load perl/intel/5.24.0

val=$SLURM_ARRAY_TASK_ID
file=`sed -n ${val}p toNorm.txt`

perl /scratch/cgsb/ercan/scripts/chip/COUNTwig_normalize_bytotalcounts-megahits_forfolder_median.pl $file
exit 0;
