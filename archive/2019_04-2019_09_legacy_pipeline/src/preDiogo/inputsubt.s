#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=subtractInput
#SBATCH --output=/scratch/las821/reports/slurm_inputsubt_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_inputsubt_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load perl/intel/5.24.0
module load r/intel/3.4.2

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p files.txt`
all=(`echo $params | sed 's/\s/ /g'`)

perl /scratch/cgsb/ercan/scripts/chip/COUNTwig_subtract_SE_folder.pl ${all[0]} ${all[1]} ${all[2]}

exit 0;
