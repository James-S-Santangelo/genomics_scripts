#!/bin/bash
#
# Script to run freebayes and generate vcf output 

## INPUTS AND OUTPUTS
export INPATH=../../data/raw/bam/bwa_mapping/RGmod

export OUTPATH=../../analysis/vcf/whole_genome
export OUTFILE="${OUTPATH}/wholeGenome_120samples"

export TMPDIR=./tmp
export REGION_FILE=../../resources/freebayes_wholeGenome.regions
export BAM_LIST=../../resources/freebayes_bams.list
export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5/GCA_005869975.1_AgR_To_v5_genomic.fna
export LOG_PATH=./log_files/freebayes
export LOG_FILE="${LOG_PATH}/freebayes.log"

# Create directories, if necessary
mkdir -p ${LOG_PATH} ${OUTPATH}

# Create list of BAM files from all BAMS in $INPATH
find ${INPATH} -name "*.bam" > ${BAM_LIST} 

## RUN FREEBAYES

{
echo "Running freebayes on all samples" > ${LOG_FILE}
~/Downloads/freebayes/scripts/freebayes-parallel \
    ${REGION_FILE} 40 \
    --fasta-reference ${PATH_TO_REFERENCE} \
    --bam-list ${BAM_LIST} \
    --use-best-n-alleles 4 \
    --report-monomorphic \
    --max-complex-gap 1 \
    --haplotype-length 1 \
    --genotype-qualities > "${OUTFILE}.vcf" &&

echo "" >> ${LOG_FILE}
echo "Done calling variant with freebayes. Sorting VCF" >> ${LOG_FILE}
bcftools sort "${OUTFILE}.vcf" > "${OUTFILE}_sorted.vcf" &&

echo "" >> ${LOG_FILE}
echo "Done sorting. Bgzipping and tabix indexing." >> ${LOG_FILE}
bgzip --force --threads 36 "${OUTFILE}_sorted.vcf" && tabix --force "${OUTFILE}_sorted.vcf.gz"
} 2>> ${LOG_FILE} 
