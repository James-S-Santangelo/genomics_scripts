#!/bin/bash
#
# Script to split freebayes VCF into separate invariant, INDEL, and SNP VCFs 

#export VCF_IN=../../analysis/vcf/whole_genome/wholeGenome_120samples_sorted.vcf.gz
export VCF_IN=../../analysis/vcf/whole_genome/test_sorted.vcf.gz
export SITES=(snps indels invariant mnps other)
export LOG_PATH=./log_files/split_variants
export FILE_BASE=$(basename ${VCF_IN%.vcf.gz})
export OUTPATH=../../analysis/vcf/whole_genome

function split_variant {
    
    local SITE_TYPE="${1}"
    local VCF_OUT=${OUTPATH}/${FILE_BASE}_${SITE_TYPE}.vcf.gz
    local LOG_FILE=${LOG_PATH}/split_variants_${SITE_TYPE}.log 

    #echo "SITE TYPE: ${SITE_TYPE}"
    #echo "VCF_OUT: ${VCF_OUT}"
    #echo "LOG FILE: ${LOG_FILE}"
    #echo ""

    echo "Splitting ${SITE_TYPE} from ${VCF_IN}" > ${LOG_FILE}
    echo "" >> ${LOG_FILE}
    {
    if [ ${SITE_TYPE} = "invariant" ]; then
        # Extract invariat sites (i.e., those with no ALT alleles)
        bcftools view --output-type z --threads 16 --include 'N_ALT = 0' ${VCF_IN} > ${VCF_OUT} 2>> ${LOG_FILE}
    elif [ ${SITE_TYPE} = "snps" ]; then
        # Extract SNPs.
        # Need to split complex events into constituent variants
        # Need to re-select SNPs following splitting
        {
        bcftools view --output-type v --types ${SITE_TYPE} ${VCF_IN} |
            vcfallelicprimitives --keep-info --keep-geno |
            bcftools view --output-type z --threads 16 --types ${SITE_TYPE} --min-alleles 2 --max-alleles 2 > ${VCF_OUT}
        } 2>> ${LOG_FILE}
    else
        # Extract other site-types (e.g., indels, mnps, other)
        bcftools view --output-type z --threads 16 --types ${SITE_TYPE} ${VCF_IN} > ${VCF_OUT} 2>> ${LOG_FILE}
    fi
} &&

echo "Tabix indexing ${VCF_OUT}" >> ${LOG_FILE}
tabix ${VCF_OUT} &&

echo "Done splitting variants and indexing VCFs!" >> ${LOG_FILE}
}

mkdir -p ${OUTPATH} ${LOG_PATH}
export -f split_variant

# Split variants in parallel across sites in SITES
parallel -j 5 split_variant ::: ${SITES[@]}  
