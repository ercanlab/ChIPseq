#!/bin/bash
#
#BATCH --verbose
#SBATCH --job-name=plotTssAlignment
#SBATCH --time=20:00
#SBATCH --nodes=1
#SBATCH --mem=10GB

module purge
module load deeptools/3.0.2

# # Get chip and input files for coverage
val=$SLURM_ARRAY_TASK_ID
params=$(sed -n ${val}p forTssAlignment.txt)
read -r reference_path score_file <<< $params

echo "reference file: $reference_path"

reference_file="$(basename $reference_path)"
score_file_prefix=${score_file%%.bw}

WORKING_DIR="${score_file_prefix}_folder"
mkdir $WORKING_DIR
mv $score_file $WORKING_DIR
cd $WORKING_DIR

chromosomes=(I II III IV V X)
printf "Computing TSS alignment profiles\n"
references=""
for chr in ${chromosomes[*]}; do
  reference="chr${chr}.bed"
  egrep "\bCHROMOSOME_${chr}\b" $reference_path | sed -r -e 's/chromosome_(I*V*|M|X)/chr\1/Ig' -e 's/MtDNA/M/gI' > $reference
  references="${references} ${reference}"
done;

references=${references# }
echo $references
computeMatrix reference-point \
  --referencePoint TSS \
  -a 2000 -b 2000\
  -R $references \
  -S $score_file \
  --skipZeros \
  -o matrix_${score_file_prefix}.gz

printf "## chr $chr: ploting profiles"
plotProfile -m matrix_${score_file_prefix}.gz\
  -out $score_file_prefix.png\
  --numPlotsPerRow 2

cd ..
mv ${WORKING_DIR}/$score_file $WORKING_DIR/$score_file_prefix.png .
rm -rf $WORKING_DIR
