#!/bin/bash
#
# Script for trimming adapters and poly-G tails from raw reads

export INPATH=/node1nfs/archive/data/trifolium/toronto_gwsd/128.120.88.251/X202SC20030702-Z01-F001/raw_data

function run_fastp {
  
  #echo ${1}
  # Get directory name and define OUTPATH
  local OUTPATH=${1}/trimmed
  local SAMPLE=$(basename $1)

  # Get filename and define new name of trimmed file
  local BASE_IN1=$(find $1 -type f -name '*_1.fq.gz' -printf "%f\n")
  local FILENAME_IN1="${BASE_IN1%.fq.gz}"
  local OUT1="${FILENAME_IN1}_trimmed.fq.gz"

  BASE_IN2=$(find $1 -type f -name '*_2.fq.gz' -printf "%f\n")
  local FILENAME_IN2="${BASE_IN2%.fq.gz}"
  local OUT2="${FILENAME_IN2}_trimmed.fq.gz"
 
  # Echo paths for debugging
  #echo ${SAMPLE}  
  #echo ${OUTPATH}
  #echo ${BASE_IN2}
  #echo ${FILENAME_IN2}
  #echo ${OUT2}

  # Create output directory if it doesn't exist
  mkdir -p ${OUTPATH}

  # Run fastp
  fastp --in1 ${1}/${BASE_IN1} \
    --in2 ${1}/${BASE_IN2} \
    --out1 ${OUTPATH}/${OUT1} \
    --out2 ${OUTPATH}/${OUT2} \
    --unpaired1 ${OUTPATH}/${SAMPLE}_unpaired.fq.gz \
    --unpaired2 ${OUTPATH}/${SAMPLE}_unpaired.fq.gz \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --overrepresentation_analysis \
    --html ${SAMPLE}_fastp.html  
  
}

export -f run_fastp

find ${INPATH} -type d -name 's_*' | parallel -j 30 --verbose run_fastp {} 
