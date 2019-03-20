#!/bin/bash


make_track_hub(){
  # extract only the jobid from the job output
  hub="$1"
  long_label="$2"
  mail="$3"

  echo "genome ce10" > genomes.txt
  echo "trackDb ce10/trackDb.txt" >> genomes.txt

  echo "hub $hub" > hub.txt
  echo "shortLabel $hub" >> hub.txt
  echo "longLabel $long_label" >> hub.txt
  echo "genomesFile genomes.txt" >> hub.txt
  echo "email $mail" >> hub.txt

  mkdir ce10

  mv *bw *bb ce10

  cd ce10

  chmod 655 *bb *bw

  ls *bw *bb > files.txt
  perl /var/trackhub_scripts/loadBWBBToUCSC.pl files.txt

  # cd back inside main hub folder
  cd ..
}

date="$(awk '!/#/ && /date:/ {print $2}' config_trackhub.yaml)"
mail="$(awk '!/#/ && /mail:/ {print $2}' config_trackhub.yaml)"
echo "date: $date"
echo "mail: $mail"
echo

hub_replicates="$date"
long_label_replicates="${date} single replicates"

replicates_dir="single_replicates"
# create hub for replicates
cd $replicates_dir

echo "Creating track hub for replicates with date ${date}..."

replicates_hub_path="/var/www/html/myHubs_v2/replicates/$date"
if [[ -d $replicates_hub_path ]]; then
  echo "Error: replicates hub with date $date already exists in $replicates_hub_path. See in the documentation how to add new replicates to this hub if that is what you want. Otherwise delete the $replicates_hub_path folder before running this command again."
  echo "Exiting..."
  exit 1;
fi

make_track_hub "$hub_replicates" "$long_label_replicates" "$mail"

# cd back into forUCSC folder
cd ..

mv $replicates_dir "$replicates_hub_path"
chmod 775 "$replicates_hub_path"

echo "Track hub created and located in $replicates_hub_path "
echo
cd averages

while read -r line || [[ -n "$line" ]]; do
  read -r chip_seqs protein <<< $line
  echo "Creating track hub for averages of replicates of the $protein protein..."
  mkdir $protein
  mv *${chip_seqs}* $protein
  cd $protein
  hub_avg="$protein"
  long_label_avg="$protein"
  avg_hub_path="/var/www/html/myHubs_v2/averages/$protein"
  if [[ -d $avg_hub_path ]]; then
    echo "Error: protein $protein hub already exists in $avg_hub_path. See in the documentation how to add files to this hub if that is what you want. Otherwise delete the $avg_hub_path folder before running this command again."
    echo "Exiting..."
    exit 1;
  fi
  make_track_hub "$hub_avg" "$long_label_avg" "$mail"
  cd ..
  mv $protein $avg_hub_path
  chmod 775 "$avg_hub_path"
  echo "Track hub created and located in $avg_hub_path"
  echo
done < forUCSC.txt

cd ~/Documents

rm -rf forUCSC

























