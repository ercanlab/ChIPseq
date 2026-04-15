# ChIP-seq pipeline - Change log

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
---

## [2.0.0] - 2025-10-14
ine structure** under´`/pipeline/` with modular organization:
### Added
- **New ChIP-seq pipeline version** under `/pipeline/` with modular organization:
  - `run_chip_SE_DG_v4.sh` — pipeline for **single-end** Illumina sequencing data (1×75 bp).
  - `run_chip_PE_DG_v2.sh` — pipeline for **paired-end** sequencing data (2×75 bp).
- Introduced environment configuration files, both defines workind directories, genome index, blacklist, and blast db:
  - `config_SE.env` — defines parameters for single-end runs.
  - `config_PE.env` — includes trimming, alignment, and quality-control parameters for paired-end runs.
- Each `.env` file now provides full reproducibility and portability across environments (HPC-compatible).
- All pipeline components (scripts + configs) reside in a single `pipeline/` directory for ease of execution and maintenance.

### Changed
- Replaced the **modular multi-script system** (separate `.s` and `.pl` components like `bamCompare.s`, `doMACS2.s`, `mergeBam.s`, etc.) with a **single executable pipeline** that sequentially performs trimming (paried-end), alignment, CPM-normalization, quality control, and peak calling.
- Simplified configuration:
  - Old version required manual listing of each sample’s metadata (sample_id, protein_name, input_file, etc.).  
  - New version only requires the path to the FASTQ directory and compliance with updated **naming rules**:

    ```
    # ChIP for paired end
    <descriptive_name>_<extract>_<sequence_ID>_input_<input_ID>_<R1/2>.fastq.gz
    Example: DPY27_N2_Emb_ext40_SE100_input_SE101_R1.fastq.gz

    # Input
    input_<strain>_<stage>_<extract>_<input_ID>_<R1/2>.fastq.gz
    Example: input_JEL1197_males_L2L3_ext726_DO98_KAPA_R2.fastq.gz

    # ChIP for single end
    <descriptive_name>_<extract>_<sequence_ID>_input_<input_ID>.fastq.gz
    Example: DPY27_N2_Emb_ext40_SE100_input_SE101.fastq.gz

    # Input
    input_<strain>_<stage>_<extract>_<input_ID>.fastq.gz
    Example: input_JEL1197_males_L2L3_ext726_DO98.fastq.gz
    ```
- Optional downstream visualization steps (`getOverlappingPeaks.s`, `plotTssAlignmentProfile.s`, `loadBWBBToUCSC.pl`, `make_trackhub.sh`) are **not integrated** and can be run manually or reintroduced in future versions.
- Documentation and README files now separated by version (`archive/` vs `pipeline/`).

### Archived
- Legacy legacy modular pipeline (2019 version) `src/` and `specs/` directories moved to `archive/2019_04-2019_09_legacy_pipeline/`.
- All scripts (e.g., `run_chipseq_v5.s`, `getOverlappingPeaks.s`) archived for reference only.
- Archived `.gitignore` and `README_legacy_2019_04-2019_09.md` for full reproducibility of historical versions.
- Added `/archive/README.md` describing the 2019 legacy version.

---

## [1.0.0] — 2019-09-01
### Added
- Initial modular ChIP-seq pipeline for single-end reads with independent scripts for alignment, normalization, and visualization.
- Manual per-sample and static configuration system using `config_v1.yaml`.
- Manual documentation (`specs/`).
- First implementation of alignment, peak calling, and visualization steps.
- Included legacy utilities such as `getOverlappingPeaks.s` and `plotTssAlignmentProfile.s`.

---

### Notes
- Version 2.0.0 marks the transition from **multi-file modular design** to a **single self-contained pipeline** emphasizing simplicity and reproducibility.  
- Archived content remains frozen to preserve provenance of prior analyses.
- Future updates may reintroduce optional visualization and trackhub modules.
