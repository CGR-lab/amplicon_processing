#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r01db22@abdn.ac.uk

module load r
module load bioconductor

basedir=$(pwd)
otu=${1:-99}

Rscript R/clusterOTUs.R $otu