#!/bin/bash
#
# Script to convert freebayes VCFs to table 

# Define Paths
#export INPATH=../../analysis/vcf/whole_genome
#export INPATH=../../analysis/vcf/chr1_only/gatk
export INPATH=../../analysis/vcf/chr1_only/freebayes
export LOG_PATH=./log_files/bcftools_filter
export OUTPATH=${INPATH}/filtered

# Define Paths to files
export LOG_FILE="${LOG_PATH}/bcftools_filter.log"
export VCF_IN="${INPATH}/tenSamplesTest_chr1_snps.vcf.gz"

# Define outfiles
export PROGRAM=$(basename ${INPATH})
export FILENAME=$(basename ${VCF_IN%.vcf.gz})
export BIALLELIC_OUT="${OUTPATH}/${PROGRAM}_${FILENAME}_biallelic.vcf.gz"
export VCF_OUT="${OUTPATH}/${PROGRAM}_${FILENAME}_biallelic_softFilter.vcf.gz"

# Make directories, if they don't exist
mkdir -p ${LOG_PATH}
mkdir -p ${OUTPATH}

#echo ${VCF_IN}
#echo ${FILENAME}
#echo #{BIALLELIC_OUT}
#echo ${VCF_OUT}

{
#Create file with only biallelic SNPs
#echo "Removing multiallelic sites from ../../${VCF_IN}" > ${LOG_FILE}
#bcftools view --min-alleles 2 \
    #--max-alleles 2 \
    #--compression-level 4 \
    #-o ${BIALLELIC_OUT} \
    #-O z ${VCF_IN} &&

echo "" >> ${LOG_FILE}
echo "Filtering biallelic VCF file ../../${BIALLELIC_OUT}" >> ${LOG_FILE}
echo "" >> ${LOG_FILE}

bcftools filter \
    -sLowQual -e'QUAL < 30' ${BIALLELIC_OUT} |
    bcftools filter \
    -o ${VCF_OUT} \
    -O z \
    --mode + \
    -sFRAC_HQ_GT -e'F_PASS(GQ < 30 | GT == "mis") > 0.2'
} 2>> ${LOG_FILE}






