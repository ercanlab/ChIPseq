#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=bedClip
#SBATCH --time=15:00
#SBATCH --nodes=1
#SBATCH --mem=10GB

####
##
##  clipBeds.s
##  This applies bedClip to bed files
##  bedClip will remove any lines that do not fall within the limits of the
##  chromosome annotation. This can occur because of extension by MACS.
##
####

module purge
module load kent/intel/20170111

val=$SLURM_ARRAY_TASK_ID
file=`sed -n ${val}p files.txt`

mv $file temp${val}.bed
bedClip temp${val}.bed /scratch/cgsb/ercan/annot/forBowtie/WS220.genome $file
rm temp${val}.bed
exit 0;
