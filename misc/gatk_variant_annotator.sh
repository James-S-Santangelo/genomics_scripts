#!/bin/bash
#
# Script to add a few annotations gatk VCF 

# Define Paths
#export INPATH=../../analysis/vcf/whole_genome
export INPATH=../../analysis/vcf/chr1_only/gatk
export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5
export LOG_PATH=./log_files/variant_annotator

# Define Paths to files
export LOG_FILE="${LOG_PATH}/gatk_variant_annotator.log"
export VCF_IN="${INPATH}/tenSamplesTest_chr1.vcf.gz"

# Define outfiles
export FILENAME=$(basename ${VCF_IN%_snps.vcf.gz})
export VCF_OUT="${INPATH}/${FILENAME}_snps_reannotated.vcf.gz"

# Make directories, if they don't exist
mkdir -p ${LOG_PATH}

#echo ${VCF_IN}
#echo ${FILENAME}
#echo ${VCF_OUT}


{
echo "Creating tables from SNPS and Indels" > ${LOG_FILE}
gatk VariantAnnotator \
    --variant ${VCF_IN} \
    --output ${VCF_OUT} \
    --reference "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
    --annotation AS_RMSMappingQuality \
    --tmp-dir ./tmp \
    --verbosity INFO &&

rm ${VCF_IN} && mv ${VCF_OUT} ${VCF_IN} && tabix -f ${VCF_IN}; 
} 2>> ${LOG_FILE}


echo "" >> ${LOG_FILE}
echo "Done reannotating VCF" >> ${LOG_FILE}
