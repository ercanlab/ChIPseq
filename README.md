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

#### Note:
NYUID is a placeholder. In commands given replace any instance of NYUID with your own NYUID:

      mkdir /scratch/$NYUID/NewChIPData    #NYUID must be replaced in this directory
      Example: mkdir /scratch/mrp420/NewChIPData   #This is correct command as I have put my own NYUID (mrp420) into the directory

#### 1. Create a new directory WD for the data (WD for working directory). This can be any name. Here we have used 'NewChIPData'.

    mkdir /scratch/$NYUID/NewChIPData
    cd /scratch/$NYUID/NewChIPData
    
#### 2. Copy the data you want to analyze into your directory. These can be bam, fastq files or both. The directory of the new data will be given by Gencore once they have finshed sequencing. Otherwise if it is an older files the bam files can be found in the ercan lab scratch, or the locations of fastqs will be in the data google spreadsheets. 

    cp /scratch/cgsb/gencore/out/Ercan/DirectoryName/FileName* /scratch/$NYUID/NewChIPData/ 

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

#### 4. Copy the analysis script to the WD

      cp /scratch/cgsb/ercan/scripts/chip/slurm/run_chipseq.sh . 
      
#### 5. Copy the configuration file to the WD and edit it to add in the details of the ChIP files

      cp /scratch/cgsb/ercan/scripts/chip/slurm/config_v1.yaml . 

After copying the configuration file (config_v1.yaml) you must edit it with the information about your data. The instructions on how to do this are in the header of the file itself. Note that the input files can be either Fastq or BAM files (simply make sure that they have the correct extension in their file name - i.e. either .fastq or .bam)

It is also possible to copy the file from google drive, edit it locally on your computer, then upload it to the HPC. The configuration file (config_v1.yaml) can be found on gdrive at:
https://drive.google.com/drive/u/1/folders/0B86yNkEPp_kmcWk5UHd5ZmVIUFk  

Once edited locally on your computer you can transfer it up to the HPCC with the following command:

      scp path/to/config_v1.yaml $NYUID@prince.hpc.nyu.edu:/scratch/$NYUID/NewChIPData/
      
#### Note:
The name of the average files that are ouutput will follow the convention:

      ${average_descriptive_name}_avg_${seq_ids}_chip
 
${average_descriptive_name} is from the configuration file. You must write this as <strong>protein target</strong>, <strong>strain</strong> and <strong>stage</strong>.   
${seq_ids} is the list of the ids that make up the average. This will be compiled by the pipeline. 


#### 6. Run the ChIP-seq pipeline. This command must be run from inside the WD.

      source run_chipseq.sh

This job will take several hours. The length of time will depend on how many files you are running.  
To check if it is still running can use the squeue command:

      squeue -u $NYUID

In the output of this command look for the row named 'chip' for current status and runtime so far. The command also gives the $JOBID of the job. These can be used to access the logs of the job (i.e. the printed messages). Logs will be stored in /scratch/$NYUID/reports/slurm_chip_$JOBID.out  
As a sanity check you can look at the top of this log file to make sure the program read the configuration file correctly. You can look at this log file while the job is running:

      less /scratch/$NYUID/reports/slurm_chip_$JOBID.out

You will also receive emails that detail the start and end of each step in the pipeline. You can keep an eye on these outputs to check when the pipeline ends. 

#### 7. Once the pipeline is finished, copy the read alignments from bowtie to archive. This command must be run from inside the WD.

      cp ReadAlignments/* /archive/s/se71/ercanlab/ReadAlignments

#### 8. Check how the job performed
      
###### A. Check your emails:
You will get a series of emails for the start and end of the different jobs that the ChIP-seq pipeline spawns. Check the emails that signal the end of a job and make sure that they have a COMPLETED status (as opposed to FAILED) and that the exit code is [0-0].  
###### B. Run the find errors command:
Some errors are not reported in the emails. And this command might catch some of these. If it finds no errors this command will output nothing. If it finds errors then it outputs the file where it found the errors and the corresponding error messages (note that the file name indicates in which part of the pipeline the error occurred).

      source /scratch/cgsb/ercan/scripts/generic/find_errors.sh       

You will see something like this if there is an error. The filename indicates that there was an error while running bowtie. Specifically in this case there was an error because the job ran out of memory:
![chip_pipeline_error_example](https://github.com/ercanlab/ChIPseq/blob/master/specs/chip_pipeline_error_example.png)

###### C. Check the log files for the job as a whole

      cat /scratch/$NYUID/reports/slurm_chip_$JOBID.out
      
###### D. Check the directory structure
Once the job had completed succesfully the WD should look as below:
![chip_pipeline_directory_example](https://github.com/ercanlab/ChIPseq/blob/master/specs/chip_pipeline_directory_example.png)



### Creating a Trackhub

#### Note: The rest of the steps should be done on ERCAN2 computer located in the computer section (let’s call it the local PC).

#### 1. Copy the forUCSC folder generated by the ChIP-seq pipeline to the local PC

      scp $NYUID@prince.hpc.nyu.edu:/scratch/NYUID/NewChIPData/forUCSC ~/Documents

#### 2. Copy the required files/scripts to make tracckhubs to your forUCSC folder

      cp /var/trackhub_scripts/config_trackhub.yaml /var/trackhub_scripts/make_trackhub.sh ~/Documents/forUCSC
      
#### 3. Edit the configuration file to match your analysis 
After copying the configuration file (config_trackhub.yaml) you must edit it with the information about your data. The instructions on how to do this are in the header of the file itself.

#### 4. Change directory to forUCSC and then run the trackhub script

      cd ~/Documents/forUCSC
      source make_trackhub.sh

Wait a few seconds

#### 5. Add the track hubs to your UCSC browser

You now have several track hubs you can add to the UCSC browser.

To add a trackhub go to https://genome.ucsc.edu . There select My data > Track hubs > My hubs
and in the URL field input the urls that are detailed below.
 
The first track hub you have is the one with all the replicates. To add this trackhub use the following url:
http://ercan2.bio.nyu.edu/myHubs_v2/replicates/date/hub.txt

Where date is the date you added in the configuration file

Second, you have the track hubs of the averages. Here you will have one track hub for each protein (or sample) you chipped. Add one by one to the browser using the url: 
/var/www/html/myHubs_v2/averages/protein/hub.txt
http://ercan2.bio.nyu.edu/myHubs_v2/averages/protein/hub.txt
Where protein is the protein(s) you had added in the ChIP-seq pipeline configuration

#### 6. Add additional files to the track hubs (in case you need to)

###### A. Put the file in the correct folder

      scp $NYUID@prince.hpc.nyu.edu:/path/to/file /var/www/html/myHubs_v2/averages/protein/ce10

###### B. Rerun the track hub configuration perl script

       cd /var/www/html/myHubs_v2/averages/protein/ce10
       ls *bw *bb > files.txt
       perl loadBWBBToUCSC.pl files.txt

###### C. Add again the trackhub to the UCSC browser
Follow step 5 and using the url: /var/www/html/myHubs_v2/averages/protein/hub.txt

#### 7. Add the track hub info to the TrackHubs tab of the Ercan_Lab_Data_WS220 document in google drive
Keep the trackhub googlesheets up to date so that others can find the data that has been processed already.The link to the document is https://docs.google.com/spreadsheets/d/1xF8nNs5dqMsMv8Ot29Hhm4Yazjs5UTwtAditn3drjCc/edit#gid=21







