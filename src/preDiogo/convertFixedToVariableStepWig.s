#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=fixedToVariableStepWig
#SBATCH --output=/scratch/las821/reports/slurm_fixedToVarWig_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_fixedToVarWig_%j.err
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load r/intel/3.4.2
module load perl/intel/5.24.0

val=$SLURM_ARRAY_TASK_ID
params=`sed -n ${val}p files.txt`
args=(`echo $params | sed 's/\s/ /g'`)

#run the R script convertFixedToVariableStepWig.R to convert a sinlge bp step fixedStep wig to a variableStep wig
fixedWig="${args[0]}"
name="${fixedWig%.*}"
variableWig="${name}_variableStep.wig"
chrI="${name}_chrIvariable.txt"
chrII="${name}_chrIIvariable.txt"
chrIII="${name}_chrIIIvariable.txt"
chrIV="${name}_chrIVvariable.txt"
chrMtDNA="${name}_chrMtDNAvariable.txt"
chrV="${name}_chrVvariable.txt"
chrX="${name}_chrXvariable.txt"
variable="${name}_variable.wig"

echo "convert $fixedWig to $variableWig"

Rscript --no-save /home/las821/scripts/convertFixedToVariableStepWig.R ${fixedWig} ${chrI} ${chrII} ${chrIII} ${chrIV} ${chrMtDNA} ${chrV} ${chrX}

cat $chrI $chrII $chrIII $chrIV $chrMtDNA $chrV $chrX > $variable

rm $chrI $chrII $chrIII $chrIV $chrMtDNA $chrV $chrX

#edit variableDPY27.wig so that you replace the file header and the CHR subheaders with variableStep wig headers
header="track type=wiggle_0 name='${name}' description='Extended tag pileup from MACS version 1.4.3 20120305 for every 1 bp'"
temp="${name}_temp.wig"
temp1="${name}_temp1.wig"
temp2="${name}_temp2.wig"
temp3="${name}_temp3.wig"

echo "$header" | cat - $variable > $temp
perl -pe 's/0\tfixedStep/variableStep/g' $temp > $temp1
perl -pe 's/\sstart=1\sstep=1//g' $temp1 > $temp2
perl -pe 's/\47/\42/g' $temp2 > $temp3

mv $temp3 $variableWig

rm $temp $temp1 $temp2
rm $variable

echo "Done"

exit 0;
