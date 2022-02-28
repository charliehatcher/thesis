#!/bin/bash

#SBATCH --array=1-100
#SBATCH --job-name=test_job
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:00
#SBATCH --mem=400GB

#Set -e allows you to test the script
set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

result_dir="/~/scratch/metab1" #EDIT

cd /~/ch17450/gctb_1.0_Linux

mkdir -p ${result_dir}

./gctb --bfile /~/biobank_pruned \
--pheno /~/UKB_pheno/metabs1.txt \
--mpheno $i --bayes NS --chain-length 30000 --burn-in 10000 --wind 0.2 \
--out /~/scratch/metab1/metab2 > test.log 2>&1




