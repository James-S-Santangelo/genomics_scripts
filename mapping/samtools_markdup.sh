#!/bin/bash
#
# Mark duplicate reads using samtools

# Define Paths
#export INPATH=../../data/raw/bam/bwa_mapping
export INPATH=../../data/raw/bam/remap

#export STATS_OUT=../../data/raw/duplications_stats
export STATS_OUT=../../data/raw/remap/duplications_stats
export LOG_PATH=./log_files/mark_dups

function mark_duplicates {

    # Define inputs and outputs
    local BAM_IN=${1}
    local FILENAME=$(basename ${BAM_IN%.bam})
    local SAMPLE=${FILENAME%_merged_sorted}
    local LOG_FILE="${LOG_PATH}/${SAMPLE}_markDups.log"
    local STATS_FILE="${STATS_OUT}/${FILENAME}_dupStats.txt"
    local OUTFILE="${INPATH}/${FILENAME}_dupsMarked.bam"

    # Echos for debugging
    #echo "BAM_IN: ${BAM_IN}"
    #echo "FILENAME: ${FILENAME}"
    #echo "SAMPLE: ${SAMPLE}"
    #echo "LOG_FILE: ${LOG_FILE}"
    #echo "STATS_FILE: ${STATS_FILE}"
    #echo "OUTFILE: ${OUTFILE}"
    #echo ""

    echo "Sorting, filling mate info, and marking duplicates for ../../${BAM_IN}" > ${LOG_FILE}
    samtools collate -O ${BAM_IN} tmp/${SAMPLE} | \
        samtools fixmate -m -O bam - - | \
        samtools sort -O bam -T tmp/ -o  - | \
        samtools markdup -T tmp/ -f ${STATS_FILE} - ${OUTFILE} 2>> ${LOG_FILE} &&

    echo "" >> ${LOG_FILE}
    echo "Indexing ../../${OUTFILE}" >> ${LOG_FILE} 
    samtools index ${OUTFILE} &&

    echo "" >> ${LOG_FILE}
    echo "Removing ../../${BAM_IN} and its index file" >> ${LOG_FILE}
    rm ${BAM_IN}
    rm "${BAM_IN}.bai"

}

mkdir -p ${LOG_PATH}
mkdir -p ${STATS_OUT}
export -f mark_duplicates

find ${INPATH} -maxdepth 1 -type f -name '*.bam' | parallel -j 30 --progress --verbose  mark_duplicates {}

