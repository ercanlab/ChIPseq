#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=getBWBB
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --mail-type=ALL

module purge
module load kent/328

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p files.txt`
read file <<< "$params"

echo "Clipping the 5th column score"
awk '{if ($5 > 1000){ $5 = 1000;}; print $0;}' $file > ${file}_tmp
mv ${file}_tmp $file

echo "Converting file $file (bed6+4) to bigBed"
perl /scratch/cgsb/ercan/scripts/trackhubs/getBWBB_v2.pl $file
exit 0;
