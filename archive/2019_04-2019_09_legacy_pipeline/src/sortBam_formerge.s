#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=sortBam
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --mem=20GB

##
## Sorts and indexes a bam file using samtools sort followed by samtools index.
##
## Input:
##  - a file named forSort.txt is assumed to exist in the directory where sortBam.s was called.
##    This file must have in each line the path to a bam file to be sorted.
##
## Warning: replaces the unsorted by the sorted bam file.

module purge
module load samtools/intel/1.6

val=$SLURM_ARRAY_TASK_ID
bamfile=$(sed -n ${val}p forSort.txt)

printf "Sorting $bamfile\n"
samtools sort -o ${bamfile}_sorted $bamfile

mv ${bamfile}_sorted $bamfile

samtools index $bamfile

exit 0;
