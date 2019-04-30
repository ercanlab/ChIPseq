# ChIP-seq analysis

This repository contains all the scripts to complete ChIP-seq analysis on the prince HPC at NYU. Methodology can be generalized to other systems, but they have been tailor made to be run by members of the Ercan lab. The current pipeline is version 5.

The below instructions will outline the steps to run the ChIP-seq analysis pipeline.

### ErcanLab_ChIP-seq_analysis_v5_slurm

Author: Sevinc Ercan - se71@nyu.edu 

Implementation:  
Diogo Mesquita - dam740@nyu.edu  
Lena Street - las821@nyu.edu   
Matt Paul - matthew.paul.2006@gmail.com  

Date: 04.2019

#### Notes:
NYUID is a placeholder. In commands given replace any instance of NYUID with your own NYUID:

      mkdir /scratch/NYUID/NewChIPData    #NYUID must be replaced in this directory
      mkdir /scratch/mrp420/NewChIPData   #This is correct command as I have put my own NYUID (mrp420) into the directory

#### 1. Create a new directory WD for the data (WD for working directory)

    mkdir /scratch/NYUID/NewChIPData
    cd /scratch/NYUID/NewChIPData
    
#### 2. Copy the data you want to analyze into your directory. These can be bam, fastq files or both. The directory of the new data will be given by Gencore once they have finshed sequencing. Otherwise if it is an older files the bam files can be found in the ercan lab scratch, or the locations of fastqs will be in the data google spreadsheets. 

    cp /scratch/cgsb/gencore/out/Ercan/DirectoryName/FileName* /scratch/NYUID/NewChIPData/ 

#### 3. Rename the files to fit with the ercan lab naming system. Ensure there are NO SPACES (if you are using BAM files instead of Fastq files, in the following examples replace .fastq by .bam)

    mv Oldfilename.fastq Newfilename.fastq

Naming Convention:	
    
    ChIP: 
    <descriptive_name>_<extract>_<sequence_ID>_input_<input_ID>.fastq
    Example: DPY27_N2_Emb_ext10_CJ132_input_CJ19.fastq

    Input:
    input_<input_id>_<strain>_<stage>_<extract>.fastq
    Example: input_CJ19_N2_Emb_ext10.fastq
    
    If the experiment involved RNAi, indicate so in the Strain name
    Example: ASH-2_N2DPY27RNAi_Emb_ext358_AKM22_input_LAS95.fastq

