#!/bin/bash
#
#BATCH --verbose
#SBATCH --job-name=bamCompare
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=30GB

##
## NOTE: This script is specific to the ChIP-seq pipeline
##

module purge
module load samtools/intel/1.6
module load deeptools/3.0.2
module load r/intel/3.4.2

compute_median(){
  bdg="$1"
  fout="$2"
  ftemp=${bdg}.tmp
  cat $bdg | awk '{print $4}' > $ftemp
  Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
          -e 'cat(median(d))' < $ftemp > $fout
  rm $ftemp
}

# set up
median_coverage_dir=MedianCoverage
mkdir -p $median_coverage_dir

# Get chip and input files for coverage
val=$SLURM_ARRAY_TASK_ID
params=$(sed -n ${val}p forBamCompare.txt)
read -r chip chip_id input input_id <<< $params

printf "\nComputing $chip coverage...\n"
chip_raw_coverage="${chip}.chip.SeqDepthNorm.bdg"
bamCoverage --bam ${chip}.bam -o $chip_raw_coverage --binSize 100 --outFileFormat bedgraph

compute_median $chip_raw_coverage $median_coverage_dir/$chip.txt
read median_chip < $median_coverage_dir/${chip}.txt
printf "\n$chip_id median coverage: $median_chip\n"

printf "\nComputing $input coverage...\n"
input_raw_coverage="${input}.input.SeqDepthNorm.bdg"
bamCoverage --bam ${input}.bam -o $input_raw_coverage --binSize 100 --outFileFormat bedgraph

compute_median $input_raw_coverage $median_coverage_dir/$input.txt
read median_input < $median_coverage_dir/${input}.txt
printf "\n$input_id median coverage: $median_input\n"

# compute inputsub
printf "\nRunning bamCompare inputsubt for $chip and $input\n"
bamCompare -b1 ${chip}.bam -b2 ${input}.bam -o ${chip}_inputsubt.bw --scaleFactors $median_input:$median_chip --operation subtract --binSize 10 --extendReads 200

printf "\nRunning bamCompare ratio for $chip and $input\n"
bamCompare -b1 ${chip}.bam -b2 ${input}.bam -o ${chip}_ratio.bw --scaleFactors $median_input:$median_chip --operation ratio --binSize 10 --extendReads 200