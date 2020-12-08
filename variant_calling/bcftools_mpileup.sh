#!/bin/bash
#
# Script to run bcftools mpileup and generate vcf output 

# Define Paths
export INPATH=../../data/raw/bam/bwa_mapping
#export INPATH=../../data/raw/subsampled_bams/RGmod

#export OUTPATH=../../analysis/vcf/whole_genome
#export OUTPATH=../../analysis/vcf/chr1_only/freebayes/
export OUTPATH=../../analysis/vcf/100genes/bcftools_mpileup/

export REGION_FILE=../../resources/100genes_sorted.bed
export BAM_LIST=../../resources/bcftools_mpileup_bams.list
export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5
export LOG_PATH=./log_files/bcftools_mpileup
export LOG_FILE="${LOG_PATH}/bcftools_mpileup.log"

mkdir -p ${LOG_PATH}
mkdir -p ${OUTPATH}

find ${INPATH} -name "*.bam" > ${BAM_LIST} 

{
echo "Running bcftools mpileup on all samples" > ${LOG_FILE}
bcftools mpileup \
    --output-type u \
    --bam-list ${BAM_LIST} \
    --fasta-ref "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
    --regions-file ${REGION_FILE} \
    --threads 36 \
    --max-depth 8000 \
    --annotate FORMAT/AD,FORMAT/DP,INFO/AD | 
    bcftools call \
    --output-type z \
    --threads 36 \
    --multiallelic-caller \
    --variants-only \
    --output "${OUTPATH}/100randomGenes_120samples.vcf.gz" -
} 2>> ${LOG_FILE} 
