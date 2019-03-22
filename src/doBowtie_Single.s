#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=bowtie
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB

module load bowtie2/intel/2.3.2
module load samtools/intel/1.6

val=$SLURM_ARRAY_TASK_ID
file=`sed -n ${val}p files.txt`

perl /scratch/cgsb/ercan/scripts/chip/slurm/getAlignment_Bowtie_SingleEnd_slurm_v2.pl $file

bam_file=${file%%.fastq}.bam

# change chromsome naming to be able to see in UCSC
samtools view -H $bam_file | sed -r -e 's/chromosome_(I*V*|M|X)/chr\1/Ig' -e 's/MtDNA/M/gI' | samtools reheader - $bam_file > ${bam_file}_chr.bam
mv ${bam_file}_chr.bam $bam_file

exit 0;
