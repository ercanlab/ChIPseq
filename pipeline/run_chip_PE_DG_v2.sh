#!/bin/bash
#SBATCH --job-name=chip_blast
#SBATCH --cpus-per-task=26
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%j_chip_blast.out
#SBATCH --mail-user=ID@nyu.edu
#SBATCH --mail-type=END

set -euo pipefail

#======================
# Don't forget
#======================
# Add your email on SBATCH --mail

# =====================
# 1. LOAD CONFIGURATION
# =====================

CONFIG="$1"

if [[ ! -f "$CONFIG" ]]; then
  echo "[ERROR] Config file not found: $CONFIG"
  exit 1
fi

source "$CONFIG"

# =====================
# 2. SETUP VARIABLES
# =====================
today=$(date +%Y_%m_%d)

fastq_dir="$WORK_DIR/fastq"

BAM="$WORK_DIR/BAM"
raw_alignment_dir="$BAM/raw_alignment"
raw_aligned="$raw_alignment_dir/${today}-raw_aligned"

trimming_dir="$BAM/trimming"
trimmed="$trimming_dir/${today}-trimmed-wo-adaptors"
unpaired="$trimming_dir/${today}-unpaired"
log_dir="$trimming_dir/log"

alignment_dir="$BAM/trimmed_alignment"
aligned="$alignment_dir/${today}-aligned"
nonaligned="$alignment_dir/${today}-non-aligned"

coverage_dir="$WORK_DIR/coverage"
CPM_norm_dir="$coverage_dir/${today}-CPM_normalization"
bam_compare_dir="$coverage_dir/${today}-bam_compare"

qc_dir="$WORK_DIR/quality_control"
fingerprint_plots="$qc_dir/${today}-fingerprint"
blast_dir="$qc_dir/${today}-blast"
subsample_dir="$blast_dir/${today}-subsample_${SUBSAMPLE_SIZE}_aligned_reads_fasta"
blast_results_dir="$blast_dir/${today}-blast_results"
counts_read="$qc_dir/${today}-counts_read"

macs2_out_dir="$WORK_DIR/MACS2output"

# =====================
# 3. LOAD MODULES
# =====================

module purge
module load trimmomatic/0.39
module load bowtie2/2.4.2
module load samtools/intel/1.11
module load deeptools/3.5.0
module load bedtools/intel/2.29.2
module load blast+/2.13.0
module load seqtk/1.5


# =====================
# 4. DETECT SAMPLES
# =====================

samples=( $(find "$fastq_dir" -name "*_R1.fastq.gz" | sed 's/_R1\.fastq\.gz$//' | xargs -n1 basename | sort -u) )

# =====================
# 5. CREATE DIRECTORIES
# =====================

mkdir -p "$BAM" "$counts_read"
mkdir -p "$raw_alignment_dir" "$raw_aligned"
mkdir -p "$trimming_dir" "$trimmed" "$unpaired" "$log_dir"
mkdir -p "$alignment_dir" "$aligned" "$nonaligned"
mkdir -p "$coverage_dir" "$CPM_norm_dir" "$bam_compare_dir"
mkdir -p "$qc_dir" "$fingerprint_plots"
mkdir -p "$blast_dir" "$subsample_dir" "$blast_results_dir"
mkdir -p "$macs2_out_dir"

# =====================
# 6. ALIGN → SORT → (REHEADER) → DEDUP → MAPQ≥20 → INDEX
# =====================
# Rationale:
# - Align PE reads with Bowtie2 and keep unaligned FASTQ for contamination checks.
# - Trim adaptors and align again.
# - Convert SAM→BAM and coordinate-sort.
# - Normalize headers (e.g., ce10 → chrI…chrM) if needed.
# - Remove duplicates with samtools markdup (-l 75 for SE 1x75, -r to remove).
# - Immediately filter MAPQ≥20; final file name remains *_sorted.bam (back-compat).
# - Indexing sorted BAM file

for sample in "${samples[@]}"; do
  echo "[INFO] Aligning: $sample"

  # 6.1
  echo "[INFO] Write SAM, and capture unaligned reads"
  bowtie2 -p $SLURM_CPUS_PER_TASK \
    -x "$INDEX_PREFIX" \
    -1 "$fastq_dir/${sample}_R1.fastq.gz" \
    -2 "$fastq_dir/${sample}_R2.fastq.gz" \
    --no-unal \
    --no-discordant --no-mixed -X $MAXINS \
    -S "$raw_aligned/${sample}.sam"
  rm "$raw_aligned/${sample}.sam"

  #6.2"
  echo "[INFO] Trimming: $sample"
  java -jar /scratch/cgsb/ercan/jar/trimmomatic-0.39.jar PE \
    -threads $SLURM_CPUS_PER_TASK -phred33 \
    "$fastq_dir/${sample}_R1.fastq.gz" "$fastq_dir/${sample}_R2.fastq.gz" \
    "$trimmed/${sample}_R1_trimmed.fastq.gz" "$unpaired/${sample}_R1_unpaired.fastq.gz" \
    "$trimmed/${sample}_R2_trimmed.fastq.gz" "$unpaired/${sample}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 CROP:$CROP MINLEN:$MINLEN &> "$log_dir/${today}_${sample}_trimmomatic.log"

  #6.3
  echo "[INFO] Aligning trimmed-$sample"
  bowtie2 -p $SLURM_CPUS_PER_TASK \
    -x "$INDEX_PREFIX" \
    -1 "$trimmed/${sample}_R1_trimmed.fastq.gz" \
    -2 "$trimmed/${sample}_R2_trimmed.fastq.gz" \
    --no-unal --no-discordant --no-mixed -X "$MAXINS" \
    --un-conc-gz "$nonaligned/${sample}_unaligned_R%.fastq.gz" \
    -S "$aligned/${sample}.sam"

  # 6.4
  echo "[INFO] Convert SAM→BAM with MAPQ≥20 FILTER"
  samtools view -@ $SLURM_CPUS_PER_TASK -b "$aligned/${sample}.sam"  > "$aligned/${sample}.bam"
  rm "$aligned/${sample}.sam"

  # 6.5
  echo "[INFO] Rename chromosomes for UCSC: ${sample}"
    header=$(samtools view -H "$aligned/${sample}.bam" | \
    sed -r \
      -e 's/chromosome_(I*V*|M|X)/chr\1/Ig' \
      -e 's/MtDNA/M/gI')
  echo "$header" | samtools reheader - "$aligned/${sample}.bam" > "$aligned/${sample}_temp_ucsc.bam"
  rm "$aligned/${sample}.bam"

  # 6.6
  echo "[INFO] REMOVE PCR duplicates (name-sort → fixmate → coord-sort → markdup)"
  # Name sort
  samtools sort -n -@ 25 -m 2G -o "$aligned/${sample}_temp_namesort.bam" "$aligned/${sample}_temp_ucsc.bam"
  # fixmate Add -m: Add ms (mate score) tags. These are used by markdup to select the best reads to keep.
  samtools fixmate -m "$aligned/${sample}_temp_namesort.bam" "$aligned/${sample}_temp_fixmate.bam"
  # Coord sort
  samtools sort -@ 25 -m 2G -o "$aligned/${sample}_temp_sorted.bam" "$aligned/${sample}_temp_fixmate.bam"
  # Remove PCR duplicates
  samtools markdup -@ "${SLURM_CPUS_PER_TASK:-1}" -r "$aligned/${sample}_temp_sorted.bam" "$aligned/${sample}_temp_sorted_rmdup.bam"

  rm -f \
  "$aligned/${sample}_temp_ucsc.bam" \
  "$aligned/${sample}_temp_namesort.bam" \
  "$aligned/${sample}_temp_fixmate.bam" \
  "$aligned/${sample}_temp_sorted.bam"

  # 6.7
  echo "[INFO] Apply MAPQ≥20 filter post remove duplicates"
  samtools view -@ 25 -b -q 20 "$aligned/${sample}_temp_sorted_rmdup.bam" > "$aligned/${sample}_sorted.bam"
  rm "$aligned/${sample}_temp_sorted_rmdup.bam"

  # 6.8
  echo "[INFO] Index final BAM"
  samtools index "$aligned/${sample}_sorted.bam"

# =====================
# 7) COVERAGE TRACKS (bamCoverage, CPM)
# =====================
# Use final *_sorted.bam (already deduplicated + MAPQ≥20)

  echo "[INFO] Creating bigWig for $sample"
  bamCoverage -b "$aligned/${sample}_sorted.bam" \
    -o "$CPM_norm_dir/${sample}.bw" \
    --binSize "$BIN_SIZE" \
    --outFileFormat bigwig \
    --minMappingQuality 20 \
    --ignoreDuplicates \
    --normalizeUsing CPM \
    --exactScaling \
    --blackListFileName "$BLACKLIST" \
    --ignoreForNormalization chrM \
    -p "$SLURM_CPUS_PER_TASK"

done

# Adding warning note for empty files in aligned directory
if ! compgen -G "$aligned/*_sorted.bam" > /dev/null; then
  echo "[WARNING] No *_sorted.bam files found under: $aligned" >&2
fi

# =====================
# 8. TREATMENT VS INPUT (bamCompare + fingerprint))
# =====================
# Pair ChIP and Input by the shared `_ext<NUM>_` token in filename.
# Operate on the final *_sorted.bam files.
shopt -s nullglob
operations=("subtract" "ratio")

# 8.1 Build map: extID -> INPUT BAM (sorted)
declare -A ext_map
echo "[INFO] Building input_map (key = extNumber)"
for input_bam in "$aligned"/input_*_ext*"_sorted.bam"; do
  file=$(basename "$input_bam")
  if [[ "$file" =~ _ext([0-9]+)_ ]]; then
    extNumber="${BASH_REMATCH[1]}"
    # Store the single sorted BAM path for that extNumber
    ext_map["$extNumber"]="$input_bam"
    echo "[INFO] Loaded INPUT: $file (extNumber=$extNumber)"
  else
    echo "[WARNING] INPUT '$file' LACKS _ext<NUM>_; skipped"
  fi
done

# 8.2 Iterative over each ChIP, find matching INPUT and produce comparisons + fingerprint
echo "[INFO] Match ChIP samples with INPUT"
for chip_bam in "$aligned"/*"_sorted.bam"; do
  chip_file=$(basename "$chip_bam")

  # Skip if it’s an input BAM
  [[ "$chip_file" == input_* ]] && continue

  # Extract extNumber from ChIP filename
  if [[ "$chip_file" =~ _ext([0-9]+)_ ]]; then
    extNumber="${BASH_REMATCH[1]}"
    echo "Processing ChIP BAM: $chip_file"
    echo "  ↳ Captured extNumber: $extNumber"
  else
    echo "[WARNING] '$chip_file' lacks _ext<NUM>_: skipped."
    continue
  fi

  # Look up the input BAM by extNumber
  control_bam="${ext_map[$extNumber]:-}"
  if [[ -z "$control_bam" ]]; then
    echo "[WARNING] No Input BAM found for extNumber=$extNumber (skipping $chip_file)."
    continue
  fi

  control_file=$(basename "$control_bam")
  echo "[INFO] Matched Control BAM: $control_file"

  # Run bamCompare for each operation
  for op in "${operations[@]}"; do
    prefix="${chip_file%.bam}"
    out_bw="$bam_compare_dir/${prefix}_${op}.bw"
    echo "[INFO] bamCompare: $chip_file vs $control_file -> $(basename "$out_bw")"
    bamCompare \
      -b1 "$chip_bam"  -b2 "$control_bam" --outFileName "$out_bw" \
      --operation "$op" $( [[ $op == "ratio" ]] && echo "--pseudocount 1" ) \
      --normalizeUsing CPM --exactScaling \
      --scaleFactorsMethod None \
      --binSize "$BIN_SIZE" --minMappingQuality 20 \
      --blackListFileName "$BLACKLIST" --ignoreForNormalization chrM \
      --outFileFormat bigwig \
      --ignoreDuplicates \
      -p "${SLURM_CPUS_PER_TASK:-1}"
    echo "[INFO] Generated $out_bw"
  done

  # Run plotFingerprint for each comparison
  prefix_chip="${chip_file%.bam}"
  prefix_input="${control_file%.bam}"
  fingerprint_png="$fingerprint_plots/${prefix_chip}_vs_${prefix_input}_fingerprint.png"
  echo "[INFO] plotFingerprint: $chip_file vs $control_file"

  plotFingerprint \
    -b "$chip_bam" "$control_bam" \
    --smartLabels --minMappingQuality 20 --skipZeros \
    --numberOfSamples "$FINGERPRINT_N" \
    --binSize "$BIN_SIZE" \
    --blackListFileName "$BLACKLIST" \
    --numberOfProcessors "$SLURM_CPUS_PER_TASK" \
    --plotFile "$fingerprint_png" \
    --plotTitle "Fingerprint: ${prefix_chip} vs ${prefix_input}"

  echo "[INFO]  $fingerprint_png generated"
done

module unload deeptools
echo "[INFO] ===  ChIP vs Input comparisons FINISHED  ==="

# =====================
# 9) MACS2 PEAK CALLING (reusing INMAP from step 8)
# =====================

module load macs2/intel/2.2.7.1
echo "[INFO] === MACS2 peak calling  ==="
# Helper: extract _ext<NUM>_ token
_get_ext() {
  local f="$1"
  if [[ "$f" =~ _ext([0-9]+)_ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo ""
  fi
}

for chip_bam in "$aligned"/*"_sorted.bam"; do
  [[ -e "$chip_bam" ]] || continue
  chip_file=$(basename "$chip_bam")
  [[ "$chip_file" == input_* ]] && continue
  ext=$(_get_ext "$chip_file")
  if [[ -z "$ext" ]]; then
    echo "[WARNING] $chip_file lacks _ext<NUM>_; skipping."
    continue
  fi

  input_bam="${ext_map[$ext]:-}"  # reuse built map from step 8 (may be empty)
  if [[ -z "$input_bam" ]]; then
    echo "[WARNING] No INPUT for ext=$ext → MACS2 will run without control."
  fi

  # Call peaks (-g ce), without no model, extend by $frag  removed because Paired-end
  macs2 callpeak \
    -t "$chip_bam" ${input_bam:+-c "$input_bam"} \
    --outdir "$macs2_out_dir" \
    -f BAMPE -g ce \
    -n "${chip_file%.bam}" \
   -q 0.05

  echo "[INFO] MACS2 done for ${chip_file%.bam}" #(extsize=$frag)"
done

# =====================
# 10) CONTAMINATION CHECK (subsample unaligned → BLAST)
# =====================
for sample in "${samples[@]}"; do
  echo "[INFO] Subsampling $SUBSAMPLE_SIZE unaligned reads from $sample for BLAST"
  gunzip -c "$nonaligned/${sample}_unaligned_R1.fastq.gz" "$nonaligned/${sample}_unaligned_R2.fastq.gz" | \
    seqtk sample -s100 - $SUBSAMPLE_SIZE | seqtk seq -a - > "$subsample_dir/${sample}_subsample_${SUBSAMPLE_SIZE}.fasta"

  echo "[INFO] Running BLAST on ${sample}_subsample_${SUBSAMPLE_SIZE}.fasta"
  blastn -query "$subsample_dir/${sample}_subsample_${SUBSAMPLE_SIZE}.fasta" \
    -db "$BLAST_DB" \
    -out "$blast_results_dir/${today}_${sample}_blastn_nt.tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -num_threads $SLURM_CPUS_PER_TASK -max_target_seqs 5 -evalue 1e-10 -parse_deflines

  echo "[INFO] Extracting best hit for ${sample}_blastn_nt.tsv"
  sort --parallel=$SLURM_CPUS_PER_TASK -k1,1 -k11,11g "$blast_results_dir/${today}_${sample}_blastn_nt.tsv" | \
    awk '!seen[$1]++' > "$blast_results_dir/${today}_${sample}_blastn_best_hits.tsv"

echo "[INFO] === Contamination diagnosis completed for ${sample}"
done

# --- CLEAN TEMPORARY BLAST ARTIFACTS (directory + large TSVs) ---
rm -rf "$subsample_dir" || true
find "$blast_results_dir" -type f -name "${today}_*_blastn_nt.tsv" -delete

# =====================
# 11) FLAG COUNTS (QC TABLE)
# =====================
echo "[INFO] Generating per-flag count table"
MASTER="$counts_read/flag_counts.csv"
echo "Sample,Flag,Count,TotalReads" > "$MASTER"
found=0
for BAM in "$aligned"/*_sorted.bam; do
  [[ -e "$BAM" ]] || continue
  found=1
  SAMPLE="$(basename "$BAM" .bam)"
  TOTAL_READS="$(samtools view -c "$BAM")"
  samtools view "$BAM" | cut -f2 | sort | uniq -c | sort -nr \
    | awk -v s="$SAMPLE" -v t="$TOTAL_READS" 'BEGIN{OFS=","} {print s, $2, $1, t}' >> "$MASTER"
  echo "[INFO] Processed $BAM (TotalReads=$TOTAL_READS)"
done
if [[ $found -eq 0 ]]; then
  echo "[WARNING] No BAMs found in: $aligned/*_sorted.bam" >&2
fi
echo "[INFO] Master flag count table: $MASTER"


# =====================
# 12) FEATURE SUMMARY (BLACKLIST / DUP / ChrM)
# =====================
echo "[INFO] Summarizing read features"
OUTPUT_FILE="$counts_read/features_summary_counts.txt"
echo -e "Sample\tTotal_Reads\tTotal_Mapped\tHighQ_Reads\tBlacklist_HighQ\tBlacklist_LowQ\tDup_HighQ\tDup_LowQ\tChrM_HighQ\tChrM_LowQ" > "$OUTPUT_FILE"

for BAM in "$aligned"/*_sorted.bam; do
  [[ -e "$BAM" ]] || continue
  SAMPLE="$(basename "$BAM" .bam)"
  echo "[INFO] Processing $SAMPLE"

  TOTAL_READS=$(samtools view -c "$BAM")
  TOTAL_MAPPED=$(samtools view -c -F 4 "$BAM")
  HIGHQ_READS=$(samtools view -c -q 20 "$BAM")

  bedtools intersect -abam "$BAM" -b "$BLACKLIST" > tmp_blacklist_unsorted.bam
  samtools sort -o tmp_blacklist.bam tmp_blacklist_unsorted.bam
  samtools index tmp_blacklist.bam
  BLACK_HQ=$(samtools view -c -q 20 tmp_blacklist.bam)
  BLACK_LQ=$(samtools view tmp_blacklist.bam | awk '$5 <= 20' | wc -l)

  # After -r removal, duplicate flags are expected ~0
  DUP_TOTAL=$(samtools view -c -f 1024 "$BAM")
  DUP_HIGH=$(samtools view -c -q 20 -f 1024 "$BAM")
  DUP_LOW=$((DUP_TOTAL - DUP_HIGH))

  # Detect mitochondrial chromosome name (MT/chrM/etc.)
  MCHR=$(samtools idxstats "$BAM" | awk '$1 ~ /(^MT$|^chrM$|M)/ {print $1; exit}')
  if [[ -n "$MCHR" ]]; then
    CHR_M_HIGH=$(samtools view -c -q 20 "$BAM" "$MCHR")
    CHR_M_LOW=$(samtools view "$BAM" "$MCHR" | awk '$5 <= 20' | wc -l)
  else
    CHR_M_HIGH=NA
    CHR_M_LOW=NA
  fi
  echo -e "$SAMPLE\t$TOTAL_READS\t$TOTAL_MAPPED\t$HIGHQ_READS\t$BLACK_HQ\t$BLACK_LQ\t$DUP_HIGH\t$DUP_LOW\t$CHR_M_HIGH\t$CHR_M_LOW" >> "$OUTPUT_FILE"
done

# Cleanup temporary files
rm -f tmp_blacklist_unsorted.bam tmp_blacklist.bam tmp_blacklist.bam.bai
echo "[INFO] Feature summary written to: $OUTPUT_FILE"


# =====================
# 13. DONE
# =====================

echo "$(date): ✅ Paired_end ChIP-seq pipeline_v2 completed"
