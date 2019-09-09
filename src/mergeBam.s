#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=mergeBam
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --mem=30GB

##
## Merges a list of bam files using samtools merge.
##
## Input:
##  - a file named forMerge.txt and is assumed to exist in the directory where
##    mergeBam.s was called.
##    This file must have in each line the arguments to the "samtools merge"
##    command separated by spaces.
##    Note that you must call this job as a job array (see https://wikis.nyu.edu/display/NYUHPC/Submitting+array+jobs).
##    That is, you must pass the "--array=1-n" argument to the sbatch command,
##    where n is the number of lines of the forMerge.txt file
##
##    For example a forMerge.txt file with the following two lines:
##      foo_out.bam foo1.bam foo2.bam foo3.bam
##      bar_out.bam bar1.bam bar2.bam bar3.bam
##    Results in running "samtools merge" twice as so:
##      samtools merge foo_out.bam foo1.bam foo2.bam foo3.bam
##      samtools merge bar_out.bam bar1.bam bar2.bam bar3.bam
##

module purge
module load samtools/intel/1.6
module load bedtools/intel/2.27.1

val=$SLURM_ARRAY_TASK_ID
params=$(sed -n ${val}p forMerge.txt)
bamfiles=($params)
mergedBam=${bamfiles[0]}

printf "Merging ${bamfiles[*]:1} into $mergedBam\n"
samtools merge $mergedBam ${bamfiles[*]:1}

printf "Converting $mergedBam from BAM to BED\n"
bedtools bamtobed -i $mergedBam > ${mergedBam}.bed

exit 0;
