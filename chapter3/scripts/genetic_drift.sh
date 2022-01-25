#!/bin/bash
#PBS -t 1-100
#PBS -l nodes=1:ppn=2,walltime=200:00:00
# Change into working directory
# Execute code
#Set -e allows you to test the script
set -e

echo "Running on ${HOSTNAME}"
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd /~/001_projects/002_selection/genetic_drift

input_file="/~/001_projects/002_selection/gctb_1.0_Linux/data/output/BayesNS/ESS_SNPs_a/cpg$i.snpRes"

output_file="ESS_a/cpg$i"

Rscript run_data.R ${input_file} ${output_file} 0.15 1 

export TMPDIR=/local
