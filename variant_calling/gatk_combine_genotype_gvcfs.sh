#!/bin/bash
#
# Run haplotype caller and output a single GVCF per sample 

# Define Paths
#export INPATH=../../analysis/vcf/whole_genome
#export INPATH=../../analysis/vcf/chr1_only/gatk/gvcf
export INPATH=../../analysis/vcf/100genes/gatk/gvcf

#export OUTPATH=../../analysis/vcf/whole_genome
export OUTPATH=../../analysis/vcf/100genes/gatk/

export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5
export LOG_PATH=./log_files/combine_genotype_gvcf

mkdir -p ${LOG_PATH}
mkdir -p ${OUTPATH}

# Inputs and outputs
export GVCF_OUT="${INPATH}/all_samples_combined.g.vcf.gz"
export GVCF_LIST="../../resources/CombineGVCF_input.list"
export LOG_FILE="${LOG_PATH}/combine_genotype_gvcf.log"
export VCF_OUT="${OUTPATH}/100randomGenes_120samples.vcf.gz"

find ${INPATH} -name "*.g.vcf.gz" > ${GVCF_LIST}

{
echo "Combining GVCF files" > ${LOG_FILE}
gatk CombineGVCFs \
 --reference "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
 --variant ${GVCF_LIST} \
 --output ${GVCF_OUT} \
 --add-output-vcf-command-line true \
 --sequence-dictionary "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.dict" \
 --tmp-dir ./tmp \
 --verbosity INFO &&

echo "" >> ${LOG_FILE}
echo "Joint genotyping of combined GVCF" >> ${LOG_FILE}
gatk GenotypeGVCFs \
    --reference "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
    --variant ${GVCF_OUT} \
    --output ${VCF_OUT} \
    --standard-min-confidence-threshold-for-calling 0 \
    --max-alternate-alleles 4 \
    --add-output-vcf-command-line true \
    --sequence-dictionary "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.dict" \
    --tmp-dir ./tmp \
    --verbosity INFO; 
} 2>> ${LOG_FILE}

echo "" >> ${LOG_FILE}
echo "Done joint genotyping" >> ${LOG_FILE}



