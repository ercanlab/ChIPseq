# ChIP-seq analysis

This repository contains all the scripts to complete ChIP-seq analysis on the prince HPC at NYU. Methodology can be generalized to other systems, but they have been tailor made to be run by members of the Ercan lab. The current pipeline is version 12.1.

The below instructions will outline the steps to run the ChIP-seq analysis pipeline.

## ErcanLab_ChIP-seq_analysis_v12.1_slurm

Author: Sevinc Ercan - se71@nyu.edu

Implementation:  
Daniel Garbozo - dg4739@nyu.edu
Date: 10.2025

Yuya Zhao - yz2954@nyu.edu
Date: 02.2025
Daniel Obaji - dno214@nyu.edu
Date: ...

Diogo Mesquita - dam740@nyu.edu  
Lena Street - las821@nyu.edu   
Matt Paul - matthew.paul.2006@gmail.com  

Date: 04.2019

## ChIP pipeline flowchart overview
![chip_pipeline_overview](https://github.com/DanielGarbozo/ChIPseq/blob/master/docs/ChIP_pipeline_overview.jpg)

## Method Overview
Sequenced reads were mapped to WS220/ce10 build of the C. elegans genome using Bowtie2 version 2.3.2 (Langmead and Salzberg, 2012). Standard mapping parameters were used for mapping with Bowtie 2. Duplicates above a expected cutoff were removed using MACS2 version 2.1.1 (Zhang, et al 2008). MACS2 determined the cutoff (typically above 1 or 2 copies) using a binomial distribution test with a p-value of 1e-5. Biological replicates were merged together with samtools (Li, et al 2009). Coverage was calculated over 100bp bins for Input subtracted ChIP enrichment scores/Chip to input ratio enrichment scores using deepTools version 3.02 (Ramírez, et al 2016). Peaks were called using MACS2. Peaks were defined using lenient settings for all biological replicates separately and using the strict settings for the merge of the biological replicates. Peaks were only considered to be true peaks if they were found in the merged replicate and in the majority of the separate replicates.

* Langmead B., Salzberg S. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923
        
        
        
         
* Zhang Y., Liu T., Meyer C.A., Eeckhoute J., Johnson D.S., Bernstein B.E., Nusbaum C., Myers R.M., Brown M., Li W., Liu X.S. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137
        
        
        
              
* Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352
        
        
        
              
* Ramírez F., Ryan D.P., Grüning B., Bhardwaj V., Kilpert F., Richter A.S., Heyne S., Dündar F., Manke T. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Res. 2016 Jul 8;44(W1):W160-5. doi: 10.1093/nar/gkw257
        
        
        
        

## Tutorial

See the tutorial on here![Tutorial_link](https://github.com/DanielGarbozo/ChIPseq/blob/master/docs/Tutorial-ErcanLab_ChIP-seq_analysis_v12.1_slurm.docx.pdf)


## How to Contribute

### 1. Fork the repository
Go to the main repository on GitHub and click **“Fork”** to create your own copy under your account.

### 2. Clone your fork
```bash
git clone <SSH link of your personal fork>
cd ChIP
```

### 3. Add the main repository as upstream

- This allows you to keep your fork updated and later create pull requests.
```bash
git remote add upstream <SSH link of the main repository>
git remote -v
```

- You should see:
```bash
origin  git@github.com:DanielGarbozo/ChIPseq.git (fetch)
origin  git@github.com:DanielGarbozo/ChIPseq.git (push)
upstream        git@github.com:ercanlab/ChIPseq.git (fetch)
upstream        git@github.com:ercanlab/ChIPseq.git (push)
```
### 4. Sync your local repo with the main branch

To get the latest updates from the main repository:
```bash
git fetch upstream dev
git merge upstream/main
```
(Use ONLY dev branch for contributions)

### 5. Make your changes

Perform your edits, add new scripts, or update documentation.

Useful file operations:
```bash
git mv <old_name> <new_name>     # move or rename files
git rm <file>                    # remove files
```

Check what has changed:
```bash
git status
```

### 6. Stage and commit your changes
```bash 
git add .                      # or specify files individually
git commit -m "Update: improved ChIP-seq pipeline"
```
### 7. Push your changes to your fork

```bash
git push origin main
```
(or replace main with master depending on your branch name)

### 8. Create a Pull Request

- Go to your fork on GitHub.
- Click “Compare & pull request.”
- Add a clear title and description explaining your contribution.
- Submit it for review 

### Notes

Keep your fork synchronized regularly:
```bash
git fetch upstream
git merge upstream/main
git push origin main
```
