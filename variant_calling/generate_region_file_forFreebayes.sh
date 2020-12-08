#!/bin/bash 
#
# Script to downsample genome and generate regions of approximately equal coverage genome-wide

## INPUTS AND OUTPUTS
export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5/GCA_005869975.1_AgR_To_v5_genomic.fna
export LOG_PATH=log_files/generate_regions
export REGIONS_OUT=../../resources/freebayes_wholeGenome.regions

# One representative BAM file. This one is ~13.7X coverage genome-wide
export BAM_IN=../../data/raw/bam/bwa_mapping/RGmod/s_7_7_merged_sorted_dupsMarked_RGmod.bam

## VARIABLE PARAMETERS
export CHROM_PREFIX='^CM*'
export NUM_REGIONS=100  # per chromosome
export SAMPLE_FRAC=0.20  # SEED.FRAC


# Create array with Chromosome IDs in bam files
chrom_array=( $(samtools idxstats ${BAM_IN} | cut -f1 | grep ${CHROM_PREFIX}) )

function generate_regions {

    local CHROM=${1}
    local LOG_FILE=${LOG_PATH}/generate_regions_${CHROM}.log
    
    # Echos for debugging
    #echo "CHROMOSOME: ${CHROM}"
    #echo "LOG_FILE: ${LOG_FILE}"

    {
    samtools view -b -s ${SAMPLE_FRAC} ${BAM_IN} ${CHROM} | bamtools coverage |\
        ~/bin/coverage_to_regions.py ${PATH_TO_REFERENCE}.fai ${NUM_REGIONS} >> ${REGIONS_OUT} 
    } 2>> ${LOG_FILE}       
}

export -f generate_regions
mkdir -p ${LOG_PATH}

touch ${REGIONS_OUT}
parallel -j 16 generate_regions ::: ${chrom_array[@]}
