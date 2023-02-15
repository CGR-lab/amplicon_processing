#!/bin/bash

conda activate microbiome

function usage {
        echo "Usage: $(basename $0) [-c CUTOFF] [-f INPUT]" 2>&1
        echo '   -c   CUTOFF    percentage similarity for OTU clustering'
        echo '   -i   INPUT     input phyloseq (.Rds) object'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring=":c:i:"

while getopts ${optstring} arg; do
  case "${arg}" in
    c) cutoff="${OPTARG:-99}" ;;
    i) input="${OPTARG:-./02_out/ps.Rds}" ;;

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

Rscript R/cluster.R $cutoff $input