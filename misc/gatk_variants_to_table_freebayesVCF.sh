#!/bin/bash
#
# Script to convert freebayes VCFs to table 

# Define Paths
#export INPATH=../../analysis/vcf/whole_genome
export INPATH=../../analysis/vcf/100genes/freebayes
export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5
export LOG_PATH=./log_files/variants_to_table
export OUTPATH=../../data/raw/snp_tables

# Define Paths to files
export LOG_FILE="${LOG_PATH}/gatk_variants_to_table.log"
#export VCF_IN="${INPATH}/freebayes_100randomGenes_120samples_snps_vcftools_filtered.recode.vcf.gz"
export VCF_IN="${INPATH}/100randomGenes_120samples_snps.vcf.gz"

# Define outfiles
export FILENAME=$(basename ${VCF_IN%.vcf.gz})
export VCF_OUT="../../data/raw/snp_tables/freebayes_${FILENAME}.table"

# Make directories, if they don't exist
mkdir -p ${LOG_PATH}
mkdir -p ${OUTPATH}

#echo ${VCF_IN}
#echo ${FILENAME}
#echo ${VCF_OUT}

echo "Creating tables from SNPS and Indels" > ${LOG_FILE}
gatk VariantsToTable \
    --variant ${VCF_IN} \
    --output ${VCF_OUT} \
    --reference "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
    --show-filtered \
    -F CHROM \
    -F POS \
    -F QUAL \
    -F REF \
    -F ALT \
    -F FILTER \
    -F DP \
    -F AF \
    -F AB \
    -F SAF \
    -F SAR \
    -F SRF \
    -F SRR \
    -F MQM \
    -F MQMR \
    -F PAIRED \
    -F PAIREDR \
    -F NO-CALL \
    -F NSAMPLES \
    -F NCALLED \
    -F MULTI-ALLELIC \
    -GF GQ \
    -GF DP \
    -GF GT \
    -GF AD \
    --tmp-dir ./tmp \
    --verbosity INFO 2>> ${LOG_FILE}

echo "" >> ${LOG_FILE}
echo "Done creating tables from SNPs and VCFS" >> ${LOG_FILE}
