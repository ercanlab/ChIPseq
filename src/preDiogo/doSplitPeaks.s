#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=splitPeaks
#SBATCH --output=/scratch/las821/reports/slurm_splitPeaks_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_splitPeaks_%j.err
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load kent/328

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p files.txt`
args=(`echo $params | sed 's/\s/ /g'`)


#sort input peaks bedfile
echo "sort Bed"
bedfile="${args[0]}"
intermed="${bedfile%.*}"
sorted="${intermed}_sorted.bed"
bedSort $bedfile $sorted


# split the peaks of the bed file
echo "splitting peaks..."
wigfile="${args[1]}"
outputDir="${args[2]}"
valley="${args[3]}"
cutoff="${args[4]}"
prefix="split_v${valley}_c${cutoff}"
outputFile="${outputDir}/${prefix}.${intermed}_sorted.subpeaks.bed"
temp="${outputDir}/${prefix}_${intermed}_temp.bed"

echo "running command: java -jar /share/apps/peak_splitter/1.0/PeakSplitter.jar -v $valley -c $cutoff -x $prefix -p $sorted -w $wigfile -o $outputDir -f false"

java -jar /share/apps/peak_splitter/1.0/PeakSplitter.jar -v $valley -c $cutoff -x $prefix -p $sorted -w $wigfile -o $outputDir -f false

#edit the output file so that it is in bed format
echo "convert output to BED format"
awk '{if (NR!=1) {print}}' $outputFile > $temp
mv $temp $outputFile
sed -i 's/chr/CHR/g' $outputFile

echo "Done"

exit 0;
