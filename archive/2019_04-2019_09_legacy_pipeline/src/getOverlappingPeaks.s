#!/bin/bash
#
#BATCH --verbose
#SBATCH --job-name=overlapping_peaks
#SBATCH --time=10:00
#SBATCH --nodes=1
#SBATCH --mem=10GB

module load bedtools/intel/2.26.0

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p forOverlappingPeaks.txt`
all=(`echo $params | sed 's/\s/ /g'`)
replicates=""
for fin in ${all[*]:1}; do
  replicates="${replicates},$fin"
done
combined_file="${all[0]}"
combined_final_file="${combined_file%.bed}_final.bed"

# remove first comma
replicates="${replicates#,}"

# get final peaks
perl /scratch/cgsb/ercan/scripts/chip/getOverlappingPeaks_Replicates.pl "$combined_file" "$replicates"

# get the final summits
summit_file="${combined_file%_peaks.bed}_summits.bed"
summit_final_file="${summit_file%.bed}_final.bed"
if [[ -f  "$summit_file" ]]; then
  bedtools intersect -wa -a "$summit_file" -b "$combined_final_file" > "$summit_final_file"
fi
exit 0;
