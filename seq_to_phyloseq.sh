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

function usage {
        echo "Usage: $(basename $0) [-abcd] [-i INPUTDIR] [-f FORID] [-r REVID]" 2>&1
        echo '   -a   AOA'
        echo '   -b   AOB'
        echo '   -c   comammox'
        echo '   -d   AOA Alves'
        echo '   -e   16S EMP'
        echo '   -i   INPUTPATH   path to the sequence fastx files'
        echo '   -f   FORID      forward read ID eg. R1/_1_'
        echo '   -r   REVID      reverse read ID eg. R2/_2_'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring=":abcdei:f:r:"

while getopts ${optstring} arg; do
  case "${arg}" in
    a) amplicon="AOA" ;;
    b) amplicon="AOB" ;;
    c) amplicon="commamox" ;;
    d) amplicon="AOA_alves" ;;
    e) amplicon="16S" ;;
    i) inputdir="${OPTARG:-./01_data/00_input}" ;;
    f) fabb="${OPTARG:-R1}" ;;
    r) rabb="${OPTARG:-R2}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# Inspect OPTIND
echo "OPTIND: $OPTIND"

basedir=$(pwd)

if [ ${amplicon} = "16S" ]
then
    clip_R1=19
    clip_R2=20
else if [ ${amplicon} = "AOB" ]
    clip_R1=21
    clip_R2=22
else
    clip_R1=22
    clip_R2=21
fi

mkdir -p 01_data/01_fastqc/01_untrimmed
mkdir -p 01_data/01_fastqc/02_trimmed
mkdir -p 01_data/02_trimmed
mkdir -p 01_data/03_filtered
mkdir -p 02_out

fastqc --noextract -o 01_data/01_fastqc/01_untrimmed $inputdir/*

## For specific primers/adapter you can set the clipping manually, which will remove the adapters and clip the amplicon on the desired
## reading frame. This example includes the parameters for AOA amoA primers (23F 616R).
#for f in $inputdir/*$fabb*; do trim_galore --fastqc -a TCGTGGGCAGCGTCAGATGTGT -a2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG --clip_R1 22 --clip_R2 21 -q 15 -o 01_data/02_trimmed --paired $f ${f/$fabb/$rabb} 2>&1 | tee 01_data/02_trimmed/$f.trimming.info; mv 01_data/02_trimmed/*fastqc* 01_data/01_fastqc/02_trimmed/; mv 01_data/02_trimmed/*report* 01_data/01_fastqc/02_trimmed/; done

for f in $inputdir/*$fabb*; do trim_galore --fastqc -a X --clip_R1 ${clip_R1} --clip_R2 ${clip_R2} -q 20 -o 01_data/02_trimmed --paired $f ${f/$fabb/$rabb}; mv 01_data/02_trimmed/*fastqc* 01_data/01_fastqc/02_trimmed/; mv 01_data/02_trimmed/*report* 01_data/01_fastqc/02_trimmed/; done

Rscript R/filter.R $amplicon $fabb $rabb

Rscript R/phyloseq.R $amplicon
