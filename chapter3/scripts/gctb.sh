#!/bin/bash

#PBS -N bayesNS
#PBS -t 1-100
#PBS -l walltime=144:00:00
#PBS -l nodes=1:ppn=4
#PBS -S /bin/bash

#Set -e allows you to test the script 
set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

result_dir="~/001_projects/002_selection/gctb_1.0_Linux/data/output/BayesNS/DNAm" #EDIT

cd ~/001_projects/002_selection/gctb_1.0_Linux/

mkdir -p ${result_dir}

gctb --bfile /~/001_projects/002_selection/phenofiles/bfiles/LCT_HLA \
--extract /~/001_projects/002_selection/ldpruned/LDpruned.prune.in \
--pheno /~/001_projects/002_selection/cpg.txt \
--mpheno $i --bayes NS --chain-length 30000 --burn-in 10000 --wind 0.2 --out $result_dir/cpg$i > test.log 2>&1







