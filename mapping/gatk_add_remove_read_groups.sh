#!/bin/bash
#
# Script to add or replace read groups

# Define Paths
export INPATH=../../data/raw/bam/bwa_mapping
#export INPATH=../../data/raw/subsampled_bams
export OUTPATH=${INPATH}/RGmod
export LOG_PATH=./log_files/add_replace_read_groups
export FASTQ_PATH=/node1nfs/archive/data/trifolium/toronto_gwsd/128.120.88.251/X202SC20030702-Z01-F001/raw_data


function add_replace_read_groups {

    # Define inputs and outputs
    local BAM_IN=${1}
    local FILENAME=$(basename ${BAM_IN%.bam})
    local SAMPLE=${FILENAME%_merged_sorted*}
    local LOG_FILE="${LOG_PATH}/${SAMPLE}_add_replace_read_groups.log"
    local OUTFILE="${OUTPATH}/${FILENAME}_RGmod.bam"

    # Get Flowcell and lane info from FASTQ header
    # Uses only read 1
    # Assumes reads from 1 sample per FASTQ
    local PATH_TO_READ_1=$(find ${FASTQ_PATH}/${SAMPLE} -type f -name "*_1_trimmed.fq.gz")
    local HEADER=$(zcat ${PATH_TO_READ_1} | head -n 1)
    local FLOWCELL="$(cut -d':' -f3 <<< "${HEADER}")"
    local LANE="$(cut -d':' -f4 <<< "${HEADER}")"
    
    ## Echos for debugging
    #echo "BAM_IN: ${BAM_IN}"
    #echo "FILENAME: ${FILENAME}"
    #echo "SAMPLE: ${SAMPLE}"
    #echo "READ_1: ${PATH_TO_READ_1}"
    #echo "HEADER: ${HEADER}"
    #echo "FLOWCELL: ${FLOWCELL}"
    #echo "LANE: ${LANE}"
    #echo "OUTFILE: ${OUTFILE}"
    #echo ""

    {
    echo "" > ${LOG_FILE}
    echo "Adding/replaceing read groups for ../../${BAM_IN}" >> ${LOG_FILE}
    gatk --java-options "-Xmx30G" AddOrReplaceReadGroups \
        --INPUT ${BAM_IN} \
        --OUTPUT ${OUTFILE} \
        --RGID ${SAMPLE} \
        --RGLB LIB1 \
        --RGPL ILLUMINA \
        --RGPU ${FLOWCELL}.${LANE}.${SAMPLE} \
        --RGSM ${SAMPLE} \
        --RGCN NOVOGENE \
        --RGPM NOVASEQ.S4.A00738 \
        --TMP_DIR ./tmp &&

    echo "" >> ${LOG_FILE}
    echo "Indexing ../../${OUTFILE}" >> ${LOG_FILE}
    samtools index ${OUTFILE} 

} 2>> ${LOG_FILE}
}

mkdir -p ${LOG_PATH}
mkdir -p ${OUTPATH}
export -f add_replace_read_groups

find ${INPATH} -maxdepth 1 -type f -name '*.bam' | parallel -j 10 --progress --verbose add_replace_read_groups {}

