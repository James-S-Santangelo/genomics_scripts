#!/bin/bash
#
# Script to split jointVCF into SNPs and Indels

# Define Paths
#export INPATH=../../analysis/vcf/whole_genome
#export INPATH=../../analysis/vcf/100genes/gatk
export INPATH=../../analysis/vcf/100genes/freebayes_75gap
export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5
export LOG_PATH=./log_files/split_variants

# Define Paths to files
export LOG_FILE="${LOG_PATH}/split_variants.log"
export VCF_IN="${INPATH}/100randomGenes_120samples.vcf.gz"

# Define outfiles
export FILENAME=$(basename ${VCF_IN%.vcf.gz})
export SNPS_OUT="${INPATH}/${FILENAME}_snps.vcf.gz"
export INDELS_OUT="${INPATH}/${FILENAME}_indels.vcf.gz"

# Make directories, if they don't exist
mkdir -p ${LOG_PATH}

#echo ${VCF_IN}
#echo ${FILENAME}
#echo ${SNPS_OUT}
#echo ${INDELS_OUT}

{
echo "Splitting ../../${VCF_IN} into SNPS and Indels" > ${LOG_FILE}
gatk SelectVariants \
    --variant ${VCF_IN} \
    --output ${SNPS_OUT} \
    --reference "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
    --select-type-to-include SNP \
    --add-output-vcf-command-line true \
    --tmp-dir ./tmp \
    --verbosity INFO &

gatk SelectVariants \
    --variant ${VCF_IN} \
    --output ${INDELS_OUT} \
    --reference "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
    --select-type-to-include INDEL \
    --add-output-vcf-command-line true \
    --tmp-dir ./tmp \
    --verbosity INFO
} 2>> ${LOG_FILE}

echo "" >> ${LOG_FILE}
echo "Done splitting variants for ../../${VCF_IN}" >> ${LOG_FILE}
