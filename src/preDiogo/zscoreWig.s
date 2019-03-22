#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=zscore
#SBATCH --output=/scratch/las821/reports/slurm_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load perl/intel/5.24.0

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p files.txt`
all=(`echo $params | sed 's/\s/ /g'`)

perl  /scratch/cgsb/ercan/scripts/chip/COUNTwig_zscore_SE_folder_peaks.pl ${all[0]} ${all[1]}
exit 0;
