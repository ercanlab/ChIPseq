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
module load bedtools/intel/2.27.1
module load samtools/intel/1.6
module load macs2/intel/2.1.1

val=$SLURM_ARRAY_TASK_ID
bamfile=$(sed -n ${val}p forSort.txt)

printf "Removing duplicates\n"
#Remove all duplicates that appear more then expected. Can modify duplicate number easily by
#adding in the --keep-dup parameter to another number i.e. --keep-dup=5
macs2 filterdup -f BAM -i $bamfile -o ${bamfile}.bed

#Sed will remove any negatvie numberr and replace it with 0. Some alignments end up -ve due
#to errors at the begininning of read.s
sed -i 's/-[0-9][0-9]*/1/' ${bamfile}.bed

#Convert the bed file output back to a bam file
bedToBam -i ${bamfile}.bed -g /scratch/cgsb/ercan/annot/forBowtie/WS220.genome > $bamfile

printf "Sorting $bamfile\n"
samtools sort -o ${bamfile}_sorted $bamfile

mv ${bamfile}_sorted $bamfile

samtools index $bamfile

exit 0;
