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

#### ChIP pipeline flowchart overview
![chip_pipeline_overview](https://github.com/ercanlab/ChIPseq/blob/master/specs/chip_pipeline_v5.png)

#### Method Overview
Sequenced reads were mapped to WS220/ce10 build of the C. elegans genome using Bowtie2 version 2.3.2 (Langmead and Salzberg, 2012). Standard mapping parameters were used for mapping with Bowtie 2. Duplicates above a expected cutoff were removed using MACS2 version 2.1.1 (Zhang, et al 2008). MACS2 determined the cutoff (typically above 1 or 2 copies) using a binomial distribution test with a p-value of 1e-5. Biological replicates were merged together with samtools (Li, et al 2009). Coverage was calculated over 100bp bins for Input subtracted ChIP enrichment scores/Chip to input ratio enrichment scores using deepTools version 3.02 (Ramírez, et al 2016). Peaks were called using MACS2. Peaks were defined using lenient settings for all biological replicates separately and using the strict settings for the merge of the biological replicates. Peaks were only considered to be true peaks if they were found in the merged replicate and in the majority of the separate replicates.

* Langmead B., Salzberg S. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.
* Zhang Y., Liu T., Meyer C.A., Eeckhoute J., Johnson D.S., Bernstein B.E., Nusbaum C., Myers R.M., Brown M., Li W., Liu X.S. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137.
* Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352.
* Ramírez F., Ryan D.P., Grüning B., Bhardwaj V., Kilpert F., Richter A.S., Heyne S., Dündar F., Manke T. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Res. 2016 Jul 8;44(W1):W160-5. doi: 10.1093/nar/gkw257.


#### 1. Create a new directory WD for the data (WD for working directory). This can be any name. Here we have used 'NewChIPData'.

```sh
mkdir /scratch/$NYUID/NewChIPData
cd /scratch/$NYUID/NewChIPData
```

#### Note:
NYUID is a placeholder. In commands given replace any instance of NYUID with your own NYUID:
```sh
mkdir /scratch/$NYUID/NewChIPData    #NYUID must be replaced in this directory
```
Example:
```sh
mkdir /scratch/mrp420/NewChIPData   #This is correct command as I have put my own NYUID (mrp420) into the directory
```

#### 2. Copy the data you want to analyze into your directory. These can be bam, fastq files or both. The directory of the new data will be given by Gencore once they have finshed sequencing. Otherwise if it is an older files the bam files can be found in the ercan lab scratch, or the locations of fastqs will be in the data google spreadsheets.

```sh
cp /scratch/cgsb/gencore/out/Ercan/DirectoryName/FileName* /scratch/$NYUID/NewChIPData/
```

#### 3. Rename the files to fit with the ercan lab naming system. Ensure there are NO SPACES (if you are using BAM files instead of Fastq files, in the following examples replace .fastq by .bam)

```sh
mv Oldfilename.fastq Newfilename.fastq
```

Naming Convention:	

    ChIP:
    <descriptive_name>_<extract>_<Sequencing_ID_for_ChIP>_input_<Sequencing_ID_for_Input>.fastq
    Example: DPY27_N2_Emb_ext10_CJ132_input_CJ19.fastq

    Input:
    input_<Sequencing_ID_for_Input>_<strain>_<stage>_<extract>.fastq
    Example: input_CJ19_N2_Emb_ext10.fastq

    If the experiment involved RNAi, indicate so in the Strain name
    Example: ASH-2_N2DPY27RNAi_Emb_ext358_AKM22_input_LAS95.fastq

#### Note:
Most errors will come from typos and the use of incorrect library names. Double-check this section to ensure that this is correct. The same goes for entering the names into the config file (step 5)  
Though placeholder library identifiers and extracts numbers will work, it essential for proper archiving of the data that the right ones are used. As the pipeline automatically saves a lot of the output data i.e. BAMs etc, it is essenital that they have the right names so when others look for the processed files they can find them. All required infromtaion can be found on the lab google sheets.

#### 4. Copy the analysis script to the WD

```sh
cp /scratch/cgsb/ercan/scripts/chip/slurm/run_chipseq.sh .
```      

#### 5. Copy the configuration file to the WD and edit it to add in the details of the ChIP files

```sh
cp /scratch/cgsb/ercan/scripts/chip/slurm/config_v1.yaml .
```

After copying the configuration file (config_v1.yaml) you must edit it with the information about your data. The instructions on how to do this are in the header of the file itself. Note that the input files can be either Fastq or BAM files (simply make sure that they have the correct extension in their file name - i.e. either .fastq or .bam)

It is also possible to copy the file from google drive, edit it locally on your computer, then upload it to the HPC. The configuration file (config_v1.yaml) can be found on gdrive at:
https://drive.google.com/drive/u/1/folders/0B86yNkEPp_kmcWk5UHd5ZmVIUFk  

Once edited locally on your computer you can transfer it up to the HPC with the following command:
```sh
scp path/to/config_v1.yaml $NYUID@prince.hpc.nyu.edu:/scratch/$NYUID/NewChIPData/
```

#### Note:
The config file must contain only the details for the samples you want to run. If there are extra entries leaving them blank or using a pound sign (#) to mask them does not work. They will be parsed and confuse the metadata files resulting in errors. Instead just delete blank entries.

#### Note:
The name of the average files that are output will follow the convention:

      ${average_descriptive_name}_avg_${seq_ids}_chip

${average_descriptive_name} is from the configuration file. You must write this as <strong>protein target</strong>_<strong>strain</strong>_<strong>stage</strong>.   
${seq_ids} is the list of the ids that make up the average. This will be compiled by the pipeline.


#### 6. Run the ChIP-seq pipeline. This command must be run from inside the WD.

```sh
source run_chipseq.sh
```

This job will take several hours. The length of time will depend on how many files you are running.  
To check if it is still running can use the squeue command:

```sh
squeue -u $NYUID
```

In the output of this command look for the row named 'chip' for current status and runtime so far. The command also gives the $JOBID of the job. These can be used to access the logs of the job (i.e. the printed messages). Logs will be stored in /scratch/$NYUID/reports/slurm_chip_$JOBID.out  
As a sanity check you can look at the top of this log file to make sure the program read the configuration file correctly. You can look at this log file while the job is running:
```sh
less /scratch/$NYUID/reports/slurm_chip_$JOBID.out
```

You will also receive emails that detail the start and end of each step in the pipeline. You can keep an eye on these outputs to check when the pipeline ends.

#### 7. Once the pipeline is finished, copy the read alignments from bowtie to archive. This command must be run from inside the WD.
```sh
cp ReadAlignments/* /archive/s/se71/ercanlab/ReadAlignments
```
#### 8. Check how the job performed

###### A. Check your emails:
You will get a series of emails for the start and end of the different jobs that the ChIP-seq pipeline spawns. Check the emails that signal the end of a job and make sure that they have a COMPLETED status (as opposed to FAILED) and that the exit code is [0-0].  
###### B. Run the find errors command:
Some errors are not reported in the emails. And this command might catch some of these. If it finds no errors this command will output nothing. If it finds errors then it outputs the file where it found the errors and the corresponding error messages (note that the file name indicates in which part of the pipeline the error occurred).

```sh
source /scratch/cgsb/ercan/scripts/generic/find_errors.sh       
```

You will see something like this if there is an error. The filename indicates that there was an error while running bowtie. Specifically in this case there was an error because the job ran out of memory:
![chip_pipeline_error_example](https://github.com/ercanlab/ChIPseq/blob/master/specs/chip_pipeline_error_example.png)

###### C. Check the log files for the job as a whole

```sh
cat /scratch/$NYUID/reports/slurm_chip_$JOBID.out
```

###### D. Check the directory structure
Once the job had completed succesfully the WD should look as below:
![chip_pipeline_directory_example](https://github.com/ercanlab/ChIPseq/blob/master/specs/chip_pipeline_directory_example.png)

#### Note:
If there is an error and you successfully fix the problem, before you decide to rerun make sure you clear your working directory of any files associated with previous run. There should just be input files (fastq or bam), config file and ChIP-seq script.

Running this code will reset directory if you are in the WD. Check there are no conflicts between rm command and names of Fastq files:

```sh
mv Fastq/* .
rm -rf BAM* Fas* for* Ra* Re* Inp* Med* MACS* re* TSS_* m*
```
If you started with BAM files (or if the issue is not at the level of mapping and you do not want to run bowtie again) you can substitute moving the contents of Fastq directory with moving contents of BAM directory.

#### Note:
These 4 examples show how ways to check that the pipeline has worked. Sometimes errors may still slip through. Keep an eye out to make sure  you get all the desired output files.


### Creating a Trackhub

#### Note: The rest of the steps should be done on ERCAN2 computer located in the computer section (let’s call it the local PC).

#### 1. Copy the forUCSC folder generated by the ChIP-seq pipeline to the local PC

```sh
scp -r $NYUID@prince.hpc.nyu.edu:/scratch/NYUID/NewChIPData/forUCSC ~/Documents
```

#### 2. Copy the required files/scripts to make tracckhubs to your forUCSC folder

```sh
cp /var/trackhub_scripts/config_trackhub.yaml /var/trackhub_scripts/make_trackhub.sh ~/Documents/forUCSC
```      

#### 3. Edit the configuration file to match your analysis
After copying the configuration file (config_trackhub.yaml) you must edit it with the information about your data. The instructions on how to do this are in the header of the file itself.

You will also want to check to make sure that noone else has added trackhubs today. The storage of the files uses the date as an index, so can be overwritten if mutiple upload trackhubs on a single day. To check run the code below. 
You will need to modify the date to todays date using this format: MM_DD_YY. If the directory exists and the command returns a short list of files, wait until the followng day to upload. Or alternatively append the date with you netID in the yaml file i.e 05_29_19_mrp420. 

```sh
ls /var/www/html/myHubs_v2/replicates/date/
```

### 4. Change directory to forUCSC and then run the trackhub script

```sh
cd ~/Documents/forUCSC
source make_trackhub.sh
```

Wait a few seconds

#### 5. Add the track hubs to your UCSC browser

You now have several track hubs you can add to the UCSC browser.

To add a trackhub go to https://genome.ucsc.edu . There select My data > Track hubs > My hubs
and in the URL field input the urls that are detailed below.

The first track hub you have is the one with all the replicates. To add this trackhub use the following url:
http://ercan2.bio.nyu.edu/myHubs_v2/replicates/date/hub.txt

Where date is the date you added in the configuration file (or if you appended the date add this. Essentially it must be the same as what was used in Step 3)

Second, you have the track hubs of the averages. Here you will have one track hub for each protein (or sample) you chipped. Add one by one to the browser using the url:
/var/www/html/myHubs_v2/averages/protein/hub.txt
http://ercan2.bio.nyu.edu/myHubs_v2/averages/protein/hub.txt
Where protein is the protein(s) you had added in the ChIP-seq pipeline configuration

#### 6. Add additional files to the track hubs (in case you need to)

###### A. Put the file in the correct folder

```sh
scp $NYUID@prince.hpc.nyu.edu:/path/to/file /var/www/html/myHubs_v2/averages/protein/ce10
```

###### B. Change directory to the folder containing the files you want to add and run the add_to_existing_hub.sh script

```sh
cd /var/www/html/myHubs_v2/averages/protein/ce10
source /var/trackhub_scripts/add_to_existing_hub.sh
```

###### C. Add again the trackhub to the UCSC browser
Follow step 5 and using the url:
/var/www/html/myHubs_v2/averages/protein/hub.txt
http://ercan2.bio.nyu.edu/myHubs_v2/averages/protein/hub.txt

#### 7. Add the track hub info to the TrackHubs tab of the Ercan_Lab_Data_WS220 document in google drive
Keep the trackhub googlesheets up to date so that others can find the data that has been processed already.The link to the document is https://docs.google.com/spreadsheets/d/1xF8nNs5dqMsMv8Ot29Hhm4Yazjs5UTwtAditn3drjCc/edit#gid=21

#### 8. Add/Backup the hub to the Ercan newhubs2019 google drive folder
Copy the hub containing folder and its contents on to google drive to back the files up and to share with lab members.
