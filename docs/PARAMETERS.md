# ChIP-seq Pipeline — Single and Paired-End Workflows

This directory contains the updated **ChIP-seq pipelines**, redesigned for clarity, reproducibility, and compatibility with HPC systems.  
Each version runs as a self-contained script and is fully configured through a single `.env` file.

---

## Overview

| File | Purpose |
|------|----------|
| `run_chip_SE_DG_v4.sh` | Complete ChIP-seq pipeline for **single-end** Illumina reads (1×75 bp). |
| `run_chip_PE_DG_v2.sh` | Complete ChIP-seq pipeline for **paired-end** Aviti reads (2×75 bp). |
| `config_SE.env` | Environment configuration for SE analysis. |
| `config_PE.env` | Environment configuration for PE analysis. |

Each `.env` file defines all required paths and parameters.
Users need to **set the working directory** and ensure **FASTQ files follow the naming rules** described below.

---

## Configuration Parameters

###  Common Parameters

| Variable | Description | Example |
|-----------|-------------|----------|
| `WORK_DIR` | Directory that contains a subfolder named `fastq` with all `.fastq.gz` files. | `/scratch/dg4739/chip_test` |
| `INDEX_PREFIX` | Path prefix of Bowtie2 genome index. | `/scratch/cgsb/ercan/annot/forBowtie/c_elegans.WS220` |
| `BLACKLIST` | BED file with genomic regions to exclude during signal generation. | `/scratch/cgsb/ercan/annot/ce10-blacklist.bed` |
| `BLAST_DB` | BLAST nt database path used for contamination checks. | `/vast/work/public/genomics/ncbi/blast/db/nt/nt` |
| `BIN_SIZE` | Bin size (bp) used in `bamCoverage` and `bamCompare`. | `10` |
| `SUBSAMPLE_SIZE` | Number of reads subsampled for BLAST QC. | `30000` |

---

### Single-End Configuration — `config_SE.env`

| Variable | Description | Example |
|-----------|-------------|----------|
| *(inherits all common parameters)* | | |
| — | — | — |
| **Notes:** | Data type: *single-end Illumina (1×75 bp)*. <br> Duplicate removal performed using `samtools markdup -l 75 -r`. <br> BAMs are filtered with MAPQ ≥ 20 before peak calling. | |

---

### Paired-End Configuration — `config_PE.env`

| Variable | Description | Example |
|-----------|-------------|----------|
| `ADAPTERS` | Directory containing adapter sequences for Trimmomatic. | `/scratch/cgsb/ercan/trimmomatic_fasta_files` |
| `CROP` | Read cropping length (bp). | `75` |
| `MINLEN` | Minimum read length after trimming. | `35` |
| `MAXINS` | Maximum allowed insert size for Bowtie2 alignment (`-X`). | `800` |
| `FINGERPRINT_N` | Number of reads used for fingerprint QC. | `2000000` |
| *(inherits all common parameters)* | | |
| **Notes:** | Data type: *paired-end Aviti (2×75 bp)*. <br> BAMs are deduplicated and filtered (MAPQ ≥ 20). | |

---

## FASTQ Naming Rules

Each input file must follow the **standardized naming convention**:

### Paired end
    ```
    # ChIP samples
    <descriptive_name>_<extract>_<sequence_ID>_input_<input_ID>_<R1/2>.fastq.gz
    Example: DPY27_N2_Emb_ext40_SE100_input_SE101_R1.fastq.gz

    # Input samples
    input_<strain>_<stage>_<extract>_<input_ID>_<R1/2>.fastq.gz
    Example: input_JEL1197_males_L2L3_ext726_DO98_KAPA_R2.fastq.gz
    ```
### Single end

    ```
    # ChIP samples
    <descriptive_name>_<extract>_<sequence_ID>_input_<input_ID>.fastq.gz
    Example: DPY27_N2_Emb_ext40_SE100_input_SE101.fastq.gz

    # Input samples
    input_<strain>_<stage>_<extract>_<input_ID>.fastq.gz
    Example: input_JEL1197_males_L2L3_ext726_DO98.fastq.gz
    ```
Both pipelines automatically matches ChIP and Input files based on the shared `<extract>` and `<input_ID>` identifiers.

---

## Execution

Run the pipeline from the project root or any accessible path in your HPC environment.

### Single-End
```bash
sbatch run_chip_SE_DG_v4.sh config_SE.env
```
### Paired-End
```bash
sbatch run_chip_PE_DG_v4.sh config_PE.env
```
## Outputs
Each run generates the following outputs under ´$WORK_DIR´:
Each run creates date-stamped folders (´e.g., 2025_09_27-*´) for reproducibility.

### Common Structure

$WORK_DIR/  
├── fastq/ → Raw input reads  
├── alignment/ or BAM/ → Aligned BAMs, non-aligned reads, trimming logs  
├── coverage/ → CPM-normalized and ChIP/Input bigWigs  
├── MACS2output/ → Peak calls (.narrowPeak, .summits.bed)  
├── quality_control/ → Fingerprint plots, read counts, BLAST QC  
└── script/ → SLURM and runtime logs  

### Key Outputs

| Category         | Main Files                       | Description                                            |
| ---------------- | -------------------------------- | ------------------------------------------------------ |
| **Alignment**    | `*.bam`                          | Sorted, deduplicated BAMs (MAPQ ≥ 20)                  |
| **Coverage**     | `*.bw`                           | Normalized signal tracks (`bamCoverage`, `bamCompare`) |
| **Peak Calling** | `*.narrowPeak`, `*.bed`          | MACS2 peak calls and summits                           |
| **QC**           | `.png`, `.txt`, `.tsv`, `.fasta` | Fingerprint, read counts, BLAST contamination check    |
| **Logs**         | `.out`                           | Pipeline and SLURM job logs                            |

## Notes and future work
- Downstream analysis modules (e.g. getOverlappingPeaks, plotTssAlignmentProfile, or UCSC trackhub gener>
These steps can be executed manually or reintroduced in future iterations.

## Manteiner: 
Daniel Garbozo - Ercan Lab, NYU Biology (2025)



