#!/bin/bash
#
# Script to compare read mapping to the reference using both BWA and NextGenMap

# Define paths 
# Validation files will be written to same directory as bams
export INPATH=../../data/raw/bam/bwa_mapping

export OUTPATH=${INPATH}/validation_files

export LOG_PATH=./log_files/validate_bams
export TMP_DIR=./tmp

function validate_bams {

    # Define inputs and outputs
    local BAM_IN=${1}
    local FILENAME=$(basename ${BAM_IN%.bam})
    local SAMPLE=${FILENAME%_merged_sorted}
    local LOG_FILE="${LOG_PATH}/${SAMPLE}_validateBam.log"
    local OUTFILE="${OUTPATH}/${FILENAME}_validation.txt"

    # Echos for debugging
    #echo "BAM_IN: ${BAM_IN}"
    #echo "FILENAME: ${FILENAME}"
    #echo "SAMPLE: ${SAMPLE}"
    #echo "LOG_FILE: ${LOG_FILE}"
    #echo "OUTFILE: ${OUTFILE}"
    #echo ""

    echo "Validating ${BAM_IN}" > ${LOG_FILE}  
    gatk ValidateSamFile \
        --INPUT ${BAM_IN} \
        --MAX_OPEN_TEMP_FILES 64000 \
        --MODE SUMMARY \
        --OUTPUT ${OUTFILE} \
        --TMP_DIR ${TMP_DIR} \
        --VERBOSITY DEBUG 2>> ${LOG_FILE} &&

    echo "" >> ${LOG_FILE}
    echo "Done validating ${BAM_IN}" >> ${LOG_FILE}
}

export -f validate_bams

# Can only run up to 10 jobs at a time on HPCNODE due to memory consumption
find ${INPATH} -maxdepth 1 -type f -name '*.bam' | parallel -j 10 --verbose --progress validate_bams {}
