#!/bin/bash

#BATCH --verbose
#SBATCH --job-name=chip
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --mail-type=ALL


#Time
start=$(date +"%T")
echo "ChIP pipeline starts at: $start"

##
# Waits for a sbatch job to finish.
# Input: string, the message returned by the sbatch command
#        e.g: 'Submitted batch job 4424072' ID i
# Will query the queue to ask if the job ID is found. Will
# sleep progress until job is cleared from the queue.
##
wait_for_job(){
  # extract only the jobid from the job output
  jobout="$1"
  jobid="${jobout##* }"

  sleep 5
  is_running=$(squeue -u $NYUID | grep $jobid | wc -l)
  while [[ $is_running -gt 0 ]]
  do
    sleep 15
    is_running=$(squeue -u $NYUID | grep $jobid | wc -l)
  done
}

#Define lab directory for automatic output archiving
ercan_chip="/scratch/cgsb/ercan/chip"

echo "Running ChIP-seq pipeline"

#Move working directory to match that of where the pipeline was initiated
cd $WORKING_DIR

mkdir -p reports

# Parse the config file and save the metadata in the metadata_dir directory.
# This nested ifelse will get the protein names, chip IDs, input IDs etc.
metadata_dir=$WORKING_DIR/metadata_dir
forMacs=$WORKING_DIR/forMacs.txt
mkdir -p $metadata_dir
egrep -A 6 "sample_id: [0-9]+" config_v1.yaml | while read -r line ; do
  val="$(cut -d ' ' -f 2 <<< $line)"
  if [[ -n "$(grep 'sample_id:' <<< $line)" ]]; then
    val="$(cut -d ' ' -f 3 <<< $line)"
    sample_id="$val"
    printf "Sample: $sample_id\n"
  elif [[ -n "$(grep 'protein_name:' <<< $line)" ]]; then
    protein="$val"
    printf "Protein: $protein\n"
  elif [[ -n "$(grep 'average_descriptive_name:' <<< $line)" ]]; then
    average_descriptive_name="$val"
    printf "Average descriptive name: $average_descriptive_name\n"
  elif [[ -n "$(grep 'input_file:' <<< $line)" ]]; then
    input="${val%%.*}"
    printf "Input: $input\n"
  elif [[ -n "$(grep 'input_seq_id:' <<< $line)" ]]; then
    input_seq_id="$val"
    printf "Input seq id: $input_seq_id\n"
  elif [[ -n "$(grep 'chip_file:' <<< $line)" ]]; then
    chip="${val%%.*}"
    printf "Chip file: $chip\n"
  elif [[ -n "$(grep 'chip_seq_id:' <<< $line)" ]]; then
    chip_seq_id="$val"
    printf "Chip seq id: $chip_seq_id\n"

    printf ""
    printf "${chip} ${chip_seq_id} ${input} ${input_seq_id}\n" >> $metadata_dir/paired_${sample_id}.txt
    printf "${protein} ${average_descriptive_name}" > $metadata_dir/metadata_${sample_id}.txt
    printf "${chip}.bam ${input}.bam ${chip} single\n" >> $forMacs
    printf "\----------\n"
  fi
done

echo "Unzipping files..."
gunzip $WORKING_DIR/*.gz


############### Run Bowtie

#List and count many fastqs in the directory
ls *.fastq > files.txt 2> /dev/null
n=$(wc -l files.txt | awk '{print $1}')

#If there are fastqs then run bowtie
if (($n > 0)); then
  echo "Running Bowtie..."
  job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_bowtie_%j.out\
                  --error=$WORKING_DIR/reports/slurm_bowtie_%j.err\
                  --mail-type=ALL\
                  --mail-user=$MAIL\
                  --array=1-$n\
                  $SBATCH_SCRIPTS/doBowtie_Single.s)

  wait_for_job "$job_out"
  echo "Bowtie finished..."

  # Record keeping. Save the bowite mapping metadata
  mkdir -p ReadAlignments
  mv Read_Alignment*txt ReadAlignments

  # Get organized. Move the Fastq files into their own directory
  mkdir Fastq
  mv *fastq Fastq
fi

#Remove the bowtie input name file
rm files.txt

# Get organized. Move the mapped BAM files to their own directory
mkdir BAM
mv *bam BAM

cd $WORKING_DIR/BAM


############### Remove duplicates, sort and index bam files

#List and count many BAMs in the directory
ls *.bam > forSort.txt
n=$(wc -l forSort.txt | awk '$0=$1')

#Duplicate reads are removed using MACS. The bedfile output is converted back to BAM for
#subsequent steps. BAM files are then sorted and indexed using sam tools, to allow further
#processing of bam files downstream
echo "Sorting and indexing BAM files..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_sortBam_%j.out\
                --error=$WORKING_DIR/reports/slurm_sortBam_%j.err\
                --mail-type=ALL\
                --mail-user=$MAIL\
                --array=1-$n\
                $SBATCH_SCRIPTS/sortBam.s)

wait_for_job "$job_out"


############### Merge BAM files

# Create the input file for the merging script
forMerge=$WORKING_DIR/BAM/forMerge.txt
for sample_file in $metadata_dir/paired_*; do
  chips=()
  inps=()
  chip_seqs=""
  inp_seqs=""
  while read -r line || [[ -n "$line" ]]; do
    read chip chip_seq_id input input_seq_id <<< $line
    chips+=("${chip}.bam")
    inps+=("${input}.bam")
    [[ -n $chip_seqs ]] && chip_seqs="${chip_seqs}_${chip_seq_id}" || chip_seqs="$chip_seq_id"
    [[ -n $inp_seqs ]] && inp_seqs="${inp_seqs}_${input_seq_id}" || inp_seqs="$input_seq_id"
  done < $sample_file

  sample_basename=$(basename $sample_file)
  tmp="${sample_basename#paired_}"
  sample_id="${tmp%.txt}"
  read -r protein average_descriptive_name < $metadata_dir/metadata_${sample_id}.txt
  printf "${chip_seqs} ${protein}\n" >> $WORKING_DIR/forUCSC.txt

  chip_seqs="avg_${chip_seqs}"
  inp_seqs="avg_${inp_seqs}"
  merged_chip_prefix="${average_descriptive_name}_${chip_seqs}_chip"
  merged_inp_prefix="${average_descriptive_name}_${inp_seqs}_input"
  merged_chip_bam="${merged_chip_prefix}.bam"
  merged_inp_bam="${merged_inp_prefix}.bam"
  printf "$merged_chip_bam ${chips[*]}\n" >> $forMerge
  printf "$merged_inp_bam ${inps[*]}\n" >> $forMerge
  printf "$merged_chip_prefix ${chip_seqs} $merged_inp_prefix ${inp_seqs}\n" >> $sample_file
  printf "$merged_chip_bam $merged_inp_bam $merged_chip_prefix merged\n" >> $forMacs
done

# Count how many BAMs mergers required using merger input file
n=$(wc -l $forMerge | awk '$0=$1')
echo "Merging BAM files..."

#BAM files are merged if they are replicates using samtools
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_mergeBam_%j.out\
                --error=$WORKING_DIR/reports/slurm_mergeBam_%j.err\
                --mail-type=ALL\
                --mail-user=$MAIL\
                --array=1-$n\
                $SBATCH_SCRIPTS/mergeBam.s)

wait_for_job "$job_out"
rm $forMerge


# sort and index merged bams

# Define and count how many merged BAM files there are
tail -q -n 1 $metadata_dir/paired_* | sed 's/[[:space:]]/.bam\n/g' | sed '2~2d' > forSort.txt
n=$(wc -l forSort.txt | awk '$0=$1')

#Merged BAM files are sorted and indexed
echo "Sorting and indexing merged BAM files..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_sortBam_%j.out\
                --error=$WORKING_DIR/reports/slurm_sortBam_%j.err\
                --mail-type=ALL\
                --mail-user=$MAIL\
                --array=1-$n\
                $SBATCH_SCRIPTS/sortBam_formerge.s)
wait_for_job "$job_out"

rm forSort.txt

# Record keeping. Put merged BAMs in one place
cp *bam $ercan_chip/ChIPseq_bamfiles/

############## BamCompare (generate bigwig files)

# create input file for bamCompare script
forBamCompare=forBamCompare.txt
cat $metadata_dir/paired_*.txt > $forBamCompare

# Define and count how many merged BAM files there are
n=$(wc -l $forBamCompare | awk '$0=$1')

# Bam compare is used to get coverage bedgraphs of input and ChIP, along with input subtracted and ratio bigwigs
echo "Running bamCompare..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_bamCompare_%j.out\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-$n\
                 $SBATCH_SCRIPTS/bamCompare.s)

wait_for_job "$job_out"
echo "bamCompare finished..."

# clean-up to remove the input metadadta
rm $forBamCompare

# organize. output files are saved into specific directories
mkdir -p InputSubtCoverage RatioCoverage RawInputCoverage RawChipCoverage
mv *chip.SeqDepthNorm.bdg RawChipCoverage
mv *input.SeqDepthNorm.bdg RawInputCoverage
mv *_inputsubt.bw InputSubtCoverage
mv *_ratio.bw RatioCoverage
mv InputSubtCoverage RatioCoverage RawInputCoverage RawChipCoverage MedianCoverage $WORKING_DIR

cd $WORKING_DIR

# Record keeping. Create a copy of files that we want to save into lab scratch so it will be archived
cp InputSubtCoverage/* $ercan_chip/ChIPseq_INSubt
cp RatioCoverage/* $ercan_chip/ChIPseq_RatioCoverage
cp MedianCoverage/* $ercan_chip/ChIPseq_MedianCoverage

############### TSS alignment profiles

cd $WORKING_DIR/InputSubtCoverage
tail -q -n 1 $metadata_dir/paired_*.txt | awk '{print $1}' > names_merged.txt

while read -r line || [[ -n "$line" ]] ; do
   file="${line}_inputsubt.bw"
   printf "$GENE_TSS_ANNOT $file\n" >> forTssAlignment.txt
 done < names_merged.txt

n=$(wc -l forTssAlignment.txt | awk '$0=$1')

echo "Running TSS alignments..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_tss_align_%j.out\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-$n\
                 $SBATCH_SCRIPTS/plotTssAlignmentProfile.s)

wait_for_job "$job_out"
echo "Tss alignments finished..."


# clean up
rm forTssAlignment.txt names_merged.txt

# organize
mkdir -p TSS_Alignment_Profiles
mv *png TSS_Alignment_Profiles
mv TSS_Alignment_Profiles $WORKING_DIR


############### Peak calling
cd $WORKING_DIR/BAM

mv $forMacs forMacs.txt

# because we're running this job in parallel it's important to create
# this directory here to avoid errors
mkdir -p MACSoutput

#Count how many files we will calculate MACs peaks for
n=$(wc -l forMacs.txt | awk '$0=$1')

#Run MACS to identify peaks of ChIP enrichment
echo "Running MACS..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_MACS2_%j.out\
                 --error=$WORKING_DIR/reports/slurm_MACS2_%j.err\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-$n\
                 $SBATCH_SCRIPTS/doMACS2.s)

wait_for_job "$job_out"
echo "MACS finished..."

#Remove the input files
rm forMacs.txt

# organize
mv MACSoutput $WORKING_DIR/
cd $WORKING_DIR/MACSoutput

# change the extension of the peak files
for file in *.narrowPeak; do
  mv "$file" "${file%.narrowPeak}".bed
done

# make the input file for overlapping peaks
forOverlappingPeaks=forOverlappingPeaks.txt
for sample_file in $metadata_dir/paired_*.txt; do
  cat -n $sample_file | sort -nr | awk -- 'ENDFILE {printf "\n"} {printf "%s_peaks.bed ",$2}' >> $forOverlappingPeaks
done

# get overlapping peaks
n=$(wc -l $forOverlappingPeaks | awk '$0=$1')

echo "Running Overlapping peaks job..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_overlapping_peaks_%j.out\
                 --error=$WORKING_DIR/reports/slurm_overlapping_peaks_%j.err\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-$n\
                 $SBATCH_SCRIPTS/getOverlappingPeaks.s)

wait_for_job "$job_out"
echo "Overlapping peaks finished..."

# organize
rm $forOverlappingPeaks

# save peaks (MACS output + final peaks)
cp *.bed $ercan_chip/ChIPseq_MACS/

cd $WORKING_DIR
rm -rf $metadata_dir
rm BAM/*.bed

############### clip Beds

cd MACSoutput/
ls *bed > files.txt

echo "bedClip to remove lines not in the limits of the chromosome annotation..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_clipBeds_%j.out\
                --error=$WORKING_DIR/reports/slurm_clipBeds_%j.err\
                --mail-type=ALL\
                --mail-user=$MAIL\
                --array=1-$n\
                $SBATCH_SCRIPTS/clipBeds.s)

wait_for_job "$job_out"
echo "bedClip finished..."

rm files.txt

############### forUCSC
cd ../
mkdir forUCSC

cp MACSoutput/*peaks.bed MACSoutput/*peaks_final.bed InputSubtCoverage/*bw RatioCoverage/*bw forUCSC

cd forUCSC

ls *bed > files.txt
n="$(wc -l files.txt | awk '{print $1}')"

echo "Converting bed to bigbed..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_getBWBB_%j.out\
                --error=$WORKING_DIR/reports/slurm_getBWBB_%j.err\
                --mail-type=ALL\
                --mail-user=$MAIL\
                --array=1-$n\
                /scratch/cgsb/ercan/scripts/trackhubs/doGetBWBB_v2.s)

wait_for_job "$job_out"
echo "Converting bed to bb finished"

rm *.bed files.txt

mkdir averages single_replicates
mv *avg* $WORKING_DIR/forUCSC.txt averages
mv *bb *bw single_replicates

#Report time
end=$(date +"%T")
echo "ChIP pipeline finishes at: $end"
