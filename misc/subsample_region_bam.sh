#!/bin/bash
#
# Script to sample bams around regions specified in bed file

# Define Paths
#export INPATH=../../data/raw/bam/bwa_mapping/RGmod
export INPATH=../../analysis/vcf/100genes/gatk/bamout
export LOG_PATH=./log_files/subsample_region_bam
export REGION_FILE=../../resources/100genes_MAFcorr_toCheck.bed
export OUTPATH="../../data/raw/bam/100genes_subset"

function subsample_region_bam {

    # Define inputs and outputs
    local BAM_IN=${1}
    local FILENAME=$(basename ${BAM_IN%.bam})
    local SAMPLE=${FILENAME%_merged_sorted*}
    local LOG_FILE="${LOG_PATH}/${SAMPLE}_subsample_region_bam.log"
    local OUTFILE="${OUTPATH}/${SAMPLE}_subset_100genes_MAFcorr_gatkHaplotypes.bam"

    #Echos for debugging
    #echo "BAM_IN: ${BAM_IN}"
    #echo "FILENAME: ${FILENAME}"
    #echo "SAMPLE: ${SAMPLE}"
    #echo "OUTFILE: ${OUTFILE}"
    #echo ""

    {
    echo "" > ${LOG_FILE}
    echo "Subseting ${REGION} from ../../${BAM_IN}" >> ${LOG_FILE}
    samtools view -b -L ${REGION_FILE} ${BAM_IN} > ${OUTFILE} &&

    echo "" >> ${LOG_FILE}
    echo "Indexing ../../${OUTFILE}" >> ${LOG_FILE}
    samtools index ${OUTFILE} &&

    echo "Done subsetting bam" >> ${LOG_FILE}
} 2>> ${LOG_FILE}
}

mkdir -p ${LOG_PATH} ${OUTPATH}
export -f subsample_region_bam

find ${INPATH} -maxdepth 1 -type f -name '*.bam' | parallel -j 30 --progress --verbose subsample_region_bam {}

