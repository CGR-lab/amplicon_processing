#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r01db22@abdn.ac.uk

module load trimgalore
module load cutadapt
module load fastqc
module load r
module load bioconductor

basedir=$(pwd)
inputdir=${1:-01_data/00_input}
fabb=${2:-_R1}
rabb=${3:-_R2}

mkdir -p 01_data/01_fastqc/01_untrimmed
mkdir -p 01_data/01_fastqc/02_trimmed
mkdir -p 01_data/02_trimmed
mkdir -p 01_data/03_filtered
mkdir -p 02_out

fastqc --noextract -o 01_data/01_fastqc/01_untrimmed $inputdir/*

for f in $inputdir/*$fabb*; do trim_galore --fastqc -a X --clip_R1 19 --clip_R2 20 -q 20 -o 01_data/02_trimmed --paired $f ${f/$fabb/$rabb}; mv 01_data/02_trimmed/*fastqc* 01_data/01_fastqc/02_trimmed/; mv 01_data/02_trimmed/*report* 01_data/01_fastqc/02_trimmed/; done

Rscript R/filter.R $fabb $rabb

Rscript R/phyloseq.R

## Useage: sbatch 1.trim.sh <Path_to_fastx_directory> <Optional forward read filename designation> <Optional reverse read filename
## designation>
## 
## <Path_to_fastx_directory>                     The directory with all raw fastx files. Filename prefix up to the read direction indicator
##                                               will correspond to the sample name
## <Optional forward read filename designation>  Part of the filename which indicates the forward reads, defaults to "_R1_"
## <Optional reverse read filename designation>  Part of the filename which indicates the reverse reads, defaults to "_R2_"
##
## For specific primers/adapter you can set the clipping manually, which will remove the adapters and clip the amplicon on the desired
## reading frame. This example includes the parameters for AOA amoA primers.

#for f in $inputdir/*$fabb*; do trim_galore --fastqc -a TCGTGGGCAGCGTCAGATGTGT -a2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG --clip_R1 22 --clip_R2 21 -q 15 -o 01_data/02_trimmed --paired $f ${f/$fabb/$rabb} 2>&1 | tee 01_data/02_trimmed/$f.trimming.info; mv 01_data/02_trimmed/*fastqc* 01_data/01_fastqc/02_trimmed/; mv 01_data/02_trimmed/*report* 01_data/01_fastqc/02_trimmed/; done

## --clip_R1 19 --clip_R2 20 : for 16S V3/V4 EMP primers
## --clip_R1 21 --clip_R2 22 : for AOB amoA assembly
## --clip_R1 22 --clip_R2 21 : for AOA amoA gap 
## Adjust based on a reading frame check
