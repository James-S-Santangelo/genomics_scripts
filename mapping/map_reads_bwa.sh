#!/bin/bash
#
# Script to compare read mapping to the reference using BWA

# Define paths 
#export INPATH=/node1nfs/archive/data/trifolium/toronto_gwsd/128.120.88.251/X202SC20030702-Z01-F001/raw_data/
export INPATH=/node1nfs/archive/data/trifolium/toronto_gwsd/128.120.88.251/X202SC20030702-Z01-F001/test_map/

#export OUTPATH=../../data/raw/bam/bwa_mapping
export OUTPATH=../../data/raw/bam/remap

export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5/GCA_005869975.1_AgR_To_v5_genomic.fna.gz
export LOG_PATH=./log_files/bwa_mapping

function run_bwa {
    
    # Get Sample ID
    local SAMPLE=$(basename ${1})
    
    # Find read files
    local READ_1=$(find ${1} -type f -name '*_1_trimmed.fq.gz' -printf "%f\n")
    local READ_2=$(find ${1} -type f -name '*_2_trimmed.fq.gz' -printf "%f\n")
    local UNPAIRED=$(find ${1} -type f -name '*_unpaired.fq.gz' -printf "%f\n")
    
    # Full path to forward, reverse, and unpaired reads
    local PATH_TO_READ_1="${1}/trimmed/${READ_1}"
    local PATH_TO_READ_2="${1}/trimmed/${READ_2}"
    local PATH_TO_UNPAIRED="${1}/trimmed/${UNPAIRED}"
    
    # Define outfile prefix
    local OUTFILE="${OUTPATH}/${SAMPLE}"
    local LOG_FILE="${LOG_PATH}/${SAMPLE}_bwaMapping.log"

    # Get Flowcell and lane info from FASTQ header
    # Uses only read 1
    # Assumes reads from 1 sample per FASTQ
    local HEADER=$(zcat ${PATH_TO_READ_1} | head -n 1)
    local FLOWCELL="$(cut -d':' -f3 <<< "${HEADER}")"
    local LANE="$(cut -d':' -f4 <<< "${HEADER}")"

    # Echos for debugging and developement
    #echo "Input: ${1}"
    #echo "Sample: ${SAMPLE}"
    #echo "Read 1 path: ${PATH_TO_READ_1}"
    #echo "Read 2 path: ${PATH_TO_READ_2}"
    #echo "Unpaired read path: ${PATH_TO_UNPAIRED}"
    #echo "Outfile: ${OUTFILE}"
    #echo "Header: ${HEADER}"
    #echo "Flowcell: ${FLOWCELL}"
    #echo "Lane: ${LANE}"
    #echo "" 

    echo "Aligning unpaired reads for sample ${SAMPLE}" >> ${LOG_FILE}
    bwa mem ${PATH_TO_REFERENCE} ${PATH_TO_UNPAIRED} \
        -R "@RG\tID:${SAMPLE}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4.A00738\tPU:${FLOWCELL}.${LANE}.${SAMPLE}\tSM:${SAMPLE}" \
        > "${OUTFILE}_unpaired.sam" 2>> ${LOG_FILE} &&
 
    echo "" >> ${LOG_FILE}    
    echo "Converting ${OUTFILE}_unpaired.sam to sorted bamfile ${OUTFILE}_unpaired_sorted.bam" >> ${LOG_FILE}  
    samtools view -b "${OUTFILE}_unpaired.sam" | samtools sort > "${OUTFILE}_unpaired_sorted.bam" 2>> ${LOG_FILE} &&
    
    echo "" >> ${LOG_FILE}    
    echo "Aligning reads for ${READ_1} and ${READ_2} to reference genome using bwa. Storing mapping as ${OUTFILE}_paired.sam" >> ${LOG_FILE}
    bwa mem ${PATH_TO_REFERENCE} ${PATH_TO_READ_1} ${PATH_TO_READ_2} \
        -R "@RG\tID:${SAMPLE}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4.A00738\tPU:${FLOWCELL}.${LANE}.${SAMPLE}\tSM:${SAMPLE}" \
        > "${OUTFILE}_paired.sam" 2>> ${LOG_FILE} && 

    echo "" >> ${LOG_FILE}    
    echo "Converting ${OUTFILE}_paired.sam to sorted bamfile ${OUTFILE}_paired_sorted.bam" >> ${LOG_FILE}  
    samtools view -b "${OUTFILE}_paired.sam" | samtools sort > "${OUTFILE}_paired_sorted.bam" 2>> ${LOG_FILE} &&
    
    echo "" >> ${LOG_FILE}    
    echo "Merging sorted paired and unpaired bam files for ${SAMPLE}" >> ${LOG_FILE}  
    samtools merge -c -f "${OUTFILE}_merged_sorted.bam" "${OUTFILE}_unpaired_sorted.bam" "${OUTFILE}_paired_sorted.bam" 2>> ${LOG_FILE} && 
    
    echo "" >> ${LOG_FILE}
    echo "Indexing ${OUTFILE}_merged_sorted.bam" >> ${LOG_FILE}
    samtools index "${OUTFILE}_merged_sorted.bam" 2>> ${LOG_FILE} &&

    echo "" >> ${LOG_FILE}    
    echo "Removing temporary files" >> ${LOG_FILE}
    rm "${OUTFILE}_unpaired.sam"
    rm "${OUTFILE}_unpaired_sorted.bam"
    rm "${OUTFILE}_paired.sam"
    rm "${OUTFILE}_paired_sorted.bam"
}

mkdir -p ${OUTPATH}
export -f run_bwa

find ${INPATH} -type d -name 's_*' | parallel -j 30 --verbose --progress run_bwa {}
