#!/bin/bash
#
# Script to convert freebayes VCFs to table 

# Define Paths
#export INPATH=../../analysis/vcf/whole_genome
#export INPATH=../../analysis/vcf/100genes/gatk
export INPATH=../../analysis/vcf/100genes/freebayes_75gap
export LOG_PATH=./log_files/vcftools_filter
export OUTPATH=${INPATH}/vcftools_filtered

# Define Paths to files
export PROGRAM=$(basename ${INPATH})
export LOG_FILE="${LOG_PATH}/${PROGRAM}_vcftools_filter.log"
export VCF_IN="${INPATH}/100randomGenes_120samples_snps.vcf.gz"

# Define temporary files outfiles
export FILENAME=$(basename ${VCF_IN%.vcf.gz})

export BIAL_QUAL_OUT="${OUTPATH}/${PROGRAM}_${FILENAME}_biallelic_qual"
export LOW_DP_MISS="${OUTPATH}/${PROGRAM}_${FILENAME}_biallelic_qual_minDP"
export MISSING_INDV="../../resources/${PROGRAM}_missing_indv"
export HIGH_PROP_MISSING="../../resources/${PROGRAM}_high_prop_missing.ind"
export HIGH_PROP_MISS_REM="${OUTPATH}/${PROGRAM}_${FILENAME}_biallelic_qual_minDP_highPropMissRem"
export VCF_OUT="${OUTPATH}/${PROGRAM}_${FILENAME}_vcftools_filtered"

# Make directories, if they don't exist
mkdir -p ${LOG_PATH}
mkdir -p ${OUTPATH}

#echo ${VCF_IN}
#echo ${FILENAME}
#echo #{BIALLELIC_OUT}
#echo ${VCF_OUT}

{
echo "Removing multiallelic sites and sites with QUAL < 30" > ${LOG_FILE}
echo "" >> ${LOG_FILE}
vcftools --gzvcf ${VCF_IN} --minQ 30 --max-alleles 2 --recode --recode-INFO-all --out ${BIAL_QUAL_OUT} &&

echo "Recoding genotypes with low DP as missing" >> ${LOG_FILE}
echo "" >> ${LOG_FILE}
vcftools --vcf "${BIAL_QUAL_OUT}.recode.vcf" --minDP 3 --recode --recode-INFO-all --out ${LOW_DP_MISS} &&

echo "Removing individuals with high proportion of missing data" >> ${LOG_FILE} &&
echo "" >> ${LOG_FILE}
vcftools --vcf "${LOW_DP_MISS}.recode.vcf" --missing-indv --out ${MISSING_INDV} &&
mawk '$5 > 0.5' "${MISSING_INDV}.imiss" | cut -f1 > ${HIGH_PROP_MISSING} &&
vcftools --vcf "${LOW_DP_MISS}.recode.vcf" --remove ${HIGH_PROP_MISSING} --recode --recode-INFO-all --out ${HIGH_PROP_MISS_REM} &&

echo "Remove sires with >20% missing data" >> ${LOG_FILE}
echo "" >> ${LOG_FILE}
vcftools --vcf "${HIGH_PROP_MISS_REM}.recode.vcf" --max-missing 0.8 --recode --recode-INFO-all --out ${VCF_OUT} &&

echo "Done filtering ../../${VCF_IN}. Removing temperary files and bgzipping" >> ${LOG_FILE}
echo "" >> ${LOG_FILE}
bgzip "${VCF_OUT}.recode.vcf" &&

rm "${BIAL_QUAL_OUT}.recode.vcf"
rm "${LOW_DP_MISS}.recode.vcf"
rm "${HIGH_PROP_MISS_REM}.recode.vcf"
} 2>> ${LOG_FILE}






