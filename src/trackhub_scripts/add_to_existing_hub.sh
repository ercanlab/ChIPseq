#Run this script from the directory containing your ChIP bw and bb files
#ie. /var/www/html/myHubs_v2/averages/DSJ266_DPY-27/ce10/
rm trackDb.txt #removes the trackDb.txt file, otherwise uploading to UCSC will error out because of duplicate information in the file
ls *bw *bb > files.txt #Generates new files.txt file for perl script to run
perl /var/trackhub_scripts/loadBWBBToUCSC.pl files.txt #Runs script to process bw and bb files to load to UCSC
chmod 655 *bb *bw #Grants access to files
chmod 771 . #Grants access to directory
