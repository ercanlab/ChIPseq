#!/bin/bash
#
#BATCH --verbose
#SBATCH --job-name=MACS
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --mem=40GB

module load macs2/intel/2.1.1

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p forMacs.txt`
read chip input name kind <<< $params

echo "Predicting fragment size for file $chip..."
macs2 predictd -i ${chip}.bed &> ${val}.txt
fragment_size=$(grep -E -o 'predicted fragment length is [[:digit:]]+' ${val}.txt | grep -E -o '[[:digit:]]+')
echo "fragment_size: $fragment_size"

# remove unnecessary folder that the predictd command creates and helper file
rm predictd ${val}.txt

echo "Calling peaks..."

if [[ $kind == "single" ]]; then
  echo "${name} is single: q = 0.05"
  macs2 callpeak -t ${chip}.bed -c ${input}.bed --outdir MACSoutput -f BED -g ce -n $name -q 0.05 --nomodel --extsize $fragment_size
else
  echo "${name} is merged: q = 0.01"
  macs2 callpeak -t ${chip}.bed -c ${input}.bed --outdir MACSoutput -f BED -g ce -n $name -q 0.01 --nomodel --extsize $fragment_size
fi

exit 0;
