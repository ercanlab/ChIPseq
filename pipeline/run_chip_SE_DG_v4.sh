#!/bin/bash
#SBATCH --job-name=chipseq
#SBATCH --cpus-per-task=26
#SBATCH --mem=60G
#SBATCH --time=35:00:00
#SBATCH --output=slurm_%j_chipseq.out
#SBATCH --mail-user=dg4739@nyu.edu
#SBATCH --mail-type=END

set -euo pipefail

# =====================
# 1) LOAD CONFIGURATION
# =====================
# Expect a single argument: path to config_SE.env file
CONFIG="$1"
if [[ ! -f "$CONFIG" ]]; then
  echo "[ERROR] Config file not found: $CONFIG"
  exit 1
fi
# shellcheck source=/dev/null
source "$CONFIG"

# =====================
# 2) DEFINE PATHS
# =====================
today=$(date +%Y_%m_%d)

fastq_dir="$WORK_DIR/fastq"

alignment_dir="$WORK_DIR/alignment"
aligned="$alignment_dir/${today}-BAM-aligned"
nonaligned="$alignment_dir/${today}-non-aligned"

coverage_dir="$WORK_DIR/coverage"
CPM_norm_dir="$coverage_dir/${today}-CPM_normalization"
bam_compare_dir="$coverage_dir/${today}-bam_compare"

qc_dir="$WORK_DIR/quality_control"
fingerprint_plots="$qc_dir/${today}-fingerprint"
blast_dir="$qc_dir/blast"
subsample_dir="$blast_dir/${today}-subsample_${SUBSAMPLE_SIZE}_aligned_reads_fasta"
blast_results_dir="$blast_dir/${today}-blast_results"
counts_read="$qc_dir/${today}-counts_read"

macs2_out_dir="$WORK_DIR/MACS2output"

# =====================
# 3) LOAD MODULES
# =====================
module purge
module load bowtie2/2.4.2
module load samtools/intel/1.11
module load deeptools/3.5.0
module load blast+/2.13.0
module load seqtk/1.5
module load bedtools/intel/2.29.2
# module load macs2/intel/2.2.7.1  # loaded right before MACS2

# =====================
# 4) DETECT SAMPLES
# =====================
# Discover unique sample basenames from fastq/*.fastq.gz (single-end inputs)
samples=( $(find "$fastq_dir" -name "*.fastq.gz" \
  | sed 's/\.fastq\.gz$//' \
  | xargs -n1 basename \
  | sort -u) )

# =====================
# 5) CREATE OUTPUT DIRECTORIES
# =====================
mkdir -p "$alignment_dir" "$aligned" "$nonaligned"
mkdir -p "$coverage_dir" "$CPM_norm_dir" "$bam_compare_dir"
mkdir -p "$qc_dir" "$fingerprint_plots"
mkdir -p "$blast_dir" "$subsample_dir" "$blast_results_dir"
mkdir -p "$counts_read"
mkdir -p "$macs2_out_dir"

# =====================
# 6) ALIGN → SORT → (REHEADER) → DEDUP → MAPQ≥20 → INDEX
# =====================
# Rationale:
# - Align SE reads with Bowtie2 and keep unaligned FASTQ for contamination checks.
# - Convert SAM→BAM and coordinate-sort.
# - Normalize headers (e.g., ce10 → chrI…chrM) if needed.
# - Remove duplicates with samtools markdup (-l 75 for SE 1x75, -r to remove).
# - Immediately filter MAPQ≥20; final file name remains *_sorted.bam (back-compat).
for sample in "${samples[@]}"; do
  echo "[INFO] Aligning: $sample"

  # 6.1
  echo "[INFO] Write SAM, and capture unaligned reads"
  bowtie2 -p "$SLURM_CPUS_PER_TASK" \
    -x "$INDEX_PREFIX" \
    -U "$fastq_dir/${sample}.fastq.gz" \
    -S "$aligned/${sample}.sam" \
    --no-unal \
    --un-gz "$nonaligned/${sample}_unaligned.fastq.gz"

  # 6.2
  echo "[INFO] Convert SAM→BAM with MAPQ≥20 FILTER and coordinate sort"
  samtools view -@ "$SLURM_CPUS_PER_TASK" -bS -q 20 "$aligned/${sample}.sam" | \
   samtools sort -@ 25 -m 2G \
      -o "$aligned/${sample}_sorted.tmp.bam" -

  # 6.3
  echo "[INFO] Change chromosome naming to be able to see un UCSC"
  samtools reheader <(samtools view -H "$aligned/${sample}_sorted.tmp.bam" | \
     sed -r -e 's/chromosome_(I*V*|M|X)/chr\1/Ig' -e 's/MtDNA/M/gI') \
    "$aligned/${sample}_sorted.tmp.bam" > "$aligned/${sample}_sorted.named.tmp.bam"
  rm  "$aligned/${sample}_sorted.tmp.bam"

  # 6.4
  echo "[INFO] REMOVE PCR duplicates on the just-sorted BAM"
  samtools markdup -@ "${SLURM_CPUS_PER_TASK:-1}" -l 75 -r \
    "$aligned/${sample}_sorted.named.tmp.bam" \
    "$aligned/${sample}_sorted.bam"
  rm  "$aligned/${sample}_sorted.named.tmp.bam"

  # 6.5
  echo "[INFO] Index final  BAM"
  samtools index -@ "${SLURM_CPUS_PER_TASK:-1}" "$aligned/${sample}_sorted.bam"

  # 6.6 Remove the large SAM to save space
  rm  "$aligned/${sample}.sam"

  echo "[INFO] FINAL_BAM -> $aligned/${sample}_sorted.bam (MAPQ≥20 + sort + dedup + index)"
done

# =====================
# 7) COVERAGE TRACKS (bamCoverage, CPM)
# =====================
# Use final *_sorted.bam (already deduplicated + MAPQ≥20)

for BAM in "$aligned"/*_sorted.bam; do
  [[ -e "$BAM" ]] || continue
  sample="$(basename "$BAM" .bam)"
  echo "[INFO] Creating bigWig for $sample"
  bamCoverage -b "$BAM" \
    -o "$CPM_norm_dir/${sample}.bw" \
    --binSize "$BIN_SIZE" \
    --outFileFormat bigwig \
    --minMappingQuality 20 \
    --extendReads 200 \
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
# 8) TREATMENT vs INPUT (bamCompare + plotFingerprint)
# =====================
# Pair ChIP and Input by the shared `_ext<NUM>_` token in filename.
# Operate on the final *_sorted.bam files.
shopt -s nullglob
operations=("subtract" "ratio")

# 8.1 Build map: extID -> INPUT BAM (sorted)
declare -A INMAP
for input_bam in "$aligned"/input_*_ext*"_sorted.bam"; do
  f=$(basename "$input_bam")
  if [[ "$f" =~ _ext([0-9]+)_ ]]; then
    INMAP["${BASH_REMATCH[1]}"]="$input_bam"
    echo "[INFO] Loaded INPUT: $f (ext=${BASH_REMATCH[1]})"
  else
    echo "[WARNING] INPUT '$f' lacks _ext<NUM>_; skipped."
  fi
done

# 8.2 For each ChIP, find matching INPUT and produce comparisons + fingerprint
for chip_bam in "$aligned"/*"_sorted.bam"; do
  [[ -e "$chip_bam" ]] || continue
  chip_file=$(basename "$chip_bam")

  # Skip if it's an input BAM
  [[ "$chip_file" == input_* ]] && continue

  # Extract extNumber (ext) from ChIP filename
  if [[ "$chip_file" =~ _ext([0-9]+)_ ]]; then
    ext="${BASH_REMATCH[1]}"
    echo "Processing ChIP BAM: $chip_file"
    echo "  ↳ Captured extNumber: $ext"
  else
    echo "[WARNING] '$chip_file' lacks _ext<NUM>_; skipped."
    continue
  fi

  # Look up the input BAM by extNumber
  control_bam="${INMAP[$ext]:-}"
  if [[ -z "$control_bam" ]]; then
    echo "[WARNING] No INPUT(sorted) for ext=$ext → skip $chip_file."
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
      -b1 "$chip_bam" -b2 "$control_bam" --outFileName "$out_bw" \
      --operation "$op" $( [[ $op == ratio ]] && echo "--pseudocount 1" ) \
      --normalizeUsing CPM --exactScaling \
      --scaleFactorsMethod None \
      --binSize "$BIN_SIZE" --minMappingQuality 20 \
      --extendReads 200 \
      --blackListFileName "$BLACKLIST" --ignoreForNormalization chrM \
      --outFileFormat bigwig \
      --ignoreDuplicates \
      -p "${SLURM_CPUS_PER_TASK:-1}"
    echo "[INFO]   ⇢ $out_bw"
  done

  # Run plotFingerprint for each comparison
  prefix_chip="${chip_file%.bam}"
  prefix_input="${control_file%.bam}"
  fingerprint_png="$fingerprint_plots/${prefix_chip}_vs_${prefix_input}_fingerprint.png"
  echo "[INFO] plotFingerprint: $chip_file vs $control_file"
  plotFingerprint \
    -b "$chip_bam" "$control_bam" \
    --smartLabels --minMappingQuality 20 --skipZeros \
    --numberOfSamples "${FINGERPRINT_N:-500000}" \
    --binSize "$BIN_SIZE" \
    --blackListFileName "$BLACKLIST" \
    --numberOfProcessors "${SLURM_CPUS_PER_TASK:-1}" \
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

  input_bam="${INMAP[$ext]:-}"  # reuse built map from step 8 (may be empty)
  if [[ -z "$input_bam" ]]; then
    echo "[WARNING] No INPUT for ext=$ext → MACS2 will run without control."
  fi

  # Try to estimate fragment length; fallback to 200 if unavailable
  predict_log="$macs2_out_dir/${chip_file%.bam}_predictd.txt"
  frag=200
  if macs2 predictd -i "$chip_bam" &> "$predict_log"; then
    got=$(grep -Eo 'predicted fragment length is [0-9]+' "$predict_log" | grep -Eo '[0-9]+')
    [[ -n "$got" ]] && frag="$got"
  fi

  # Call peaks (-g ce), no model, extend by $frag
  macs2 callpeak \
    -t "$chip_bam" ${input_bam:+-c "$input_bam"} \
    --outdir "$macs2_out_dir" \
    -f BAM -g ce \
    -n "${chip_file%.bam}" \
    --nomodel --extsize "$frag" \
    -q 0.05

  echo "[INFO] MACS2 done for ${chip_file%.bam} (extsize=$frag)"
done

# =====================
# 10) CONTAMINATION CHECK (subsample unaligned → BLAST)
# =====================
for sample in "${samples[@]}"; do
  echo "[INFO] Subsampling $SUBSAMPLE_SIZE unaligned reads from $sample for BLAST"
  gunzip -c "$nonaligned/${sample}_unaligned.fastq.gz" | \
     seqtk sample -s100 - $SUBSAMPLE_SIZE | \
     seqtk seq -a - > "$subsample_dir/${sample}_subsample_${SUBSAMPLE_SIZE}.fasta"

  echo "[INFO] Running BLAST on ${sample}_subsample_${SUBSAMPLE_SIZE}.fasta"
  blastn -query "$subsample_dir/${sample}_subsample_${SUBSAMPLE_SIZE}.fasta" \
    -db "$BLAST_DB" \
    -out "$blast_results_dir/${today}_${sample}_blastn_nt.tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -num_threads $SLURM_CPUS_PER_TASK -max_target_seqs 5 -evalue 1e-10 -parse_deflines

  echo "[INFO] Extracting best hit for ${sample}_blastn_nt.tsv"
  sort --parallel=$SLURM_CPUS_PER_TASK -k1,1 -k11,11g "$blast_results_dir/${today}_${sample}_blastn_nt.tsv" | \
     awk '!seen[$1]++' > "$blast_results_dir/${today}_${sample}_blastn_best_hits.tsv"
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
# 13) DONE
# =====================
echo "$(date): ✅ Single-end ChIP-seq pipeline_v2 completed"
