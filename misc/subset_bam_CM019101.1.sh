#!/bin/bash
#
# Script to subset bam files for 10 test samples to include only reads
# that map to CM019101.1

# Define paths
export INPATH=../../data/raw/bam/bwa_mapping
export OUTPATH=../../data/raw/subsampled_bams

# Array with IDs of samples to subsample
samples=(s_100_13 s_116_1 s_23_12 s_37_18 s_41_13 s_43_15 s_54_9 s_7_7 s_83_10 s_97_3)

# Create array with full path to samples to subsamples
fullPathsToSamples=()
for ((i = 0; i < ${#samples[@]}; i++))
do
    fullPathsToSamples[$i]=$(find ${INPATH} -maxdepth 1 -type f -name "${samples[$i]}_*.bam")
    #echo "${fullPathsToSamples[$i]}"
done
#echo ${fullPathsToSamples[@]}

function subsample_bam {

    # Variables
    local BAM_IN=${1}
    local FILENAME=$(basename ${BAM_IN%.bam})
    local OUTFILE="${OUTPATH}/${FILENAME}_CM019101.1.bam"

    # Echos for debugging
    #echo "INPUT: ${BAM_IN}"
    #echo "FILENAME: ${FILENAME}"
    #echo "OUTFILE: ${OUTFILE}"

    echo "Subsampling ${BAM_IN}"
    samtools view -b ${BAM_IN} CM019101.1 > ${OUTFILE} &&

    echo "Indexing ${OUTFILE}"
    samtools index ${OUTFILE}
    
}

mkdir -p ${OUTPATH}
export -f subsample_bam

parallel -j 10 --ungroup subsample_bam ::: ${fullPathsToSamples[@]}
