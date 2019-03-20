#!/bin/sh

# Parse input
if [[ "$#" -gt 0 ]]; then
  echo "Error: this script takes no arguments."
  echo "Usage: ./run_chipseq.sh"
  exit -1;
fi

export WORKING_DIR=$(pwd)
export YAML_CONFIG=$WORKING_DIR/config_v1.yaml
export MAIL=$(egrep 'mail:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export NYUID=$(egrep 'nyuid:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export SBATCH_SCRIPTS=$(egrep 'sbatch_scripts:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export GENE_TSS_ANNOT=$(egrep 'gene_tss_annot:' $YAML_CONFIG | awk '!/^#/ && $0=$2')

mkdir -p /scratch/$NYUID/reports

sbatch --output=/scratch/$NYUID/reports/slurm_chip_%j.out\
       --mail-type=ALL\
       --mail-user=$MAIL\
       --export=WORKING_DIR,YAML_CONFIG,MAIL,NYUID,SBATCH_SCRIPTS,GENE_TSS_ANNOT\
       $SBATCH_SCRIPTS/run_chipseq_v5.s
