#!/bin/bash
#
# Run haplotype caller and output a single GVCF per sample 

# Define Paths
export INPATH=../../data/raw/bam/bwa_mapping/RGmod
#export INPATH=../../data/raw/subsampled_bams

#export OUTPATH=../../analysis/vcf/whole_genome
#export OUTPATH=../../analysis/vcf/chr1_only/gatk/gvcf
export OUTPATH=../../analysis/vcf/100genes/gatk/gvcf
export BAM_OUT=${OUTPATH}/../bamout

export REGIONS_FILE=../../resources/100genes_sorted.bed
export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5
export LOG_PATH=./log_files/haplotype_caller

function run_haplotypeCaller {
    
    # Define inputs and outputs
    local BAM_IN=${1}
    local FILENAME=$(basename ${BAM_IN%.bam})
    local LOG_FILE="${LOG_PATH}/${FILENAME}_gatkHC.log"
    local OUTFILE="${OUTPATH}/${FILENAME}.g.vcf.gz"

    # Echos for debugging
    #echo "BAM_IN: ${BAM_IN}"
    #echo "FILENAME: ${FILENAME}"
    #echo "LOG_FILE: ${LOG_FILE}"
    #echo "OUTFILE: ${OUTFILE}"
    #echo ""

    echo "Generating GVCF for ../../${BAM_IN}" > ${LOG_FILE}
    gatk --java-options "-Xmx30G" HaplotypeCaller \
        --input ${BAM_IN} \
        --output ${OUTFILE} \
        --reference "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.fna" \
        --intervals ${REGIONS_FILE} \
        --max-reads-per-alignment-start 0 \
        --standard-min-confidence-threshold-for-calling 0 \
        --max-alternate-alleles 4 \
        --min-base-quality-score 20 \
        --minimum-mapping-quality 30 \
        --output-mode EMIT_VARIANTS_ONLY \
        --sample-ploidy 2 \
        --add-output-vcf-command-line true \
        --tmp-dir ./tmp \
        --verbosity INFO \
        --emit-ref-confidence GVCF \
        --bam-output ${BAM_OUT}/${FILENAME}_bamout.bam \
        --sequence-dictionary "${PATH_TO_REFERENCE}/GCA_005869975.1_AgR_To_v5_genomic.dict" 2>> ${LOG_FILE}

    echo "" >> ${LOG_FILE}
    echo "Done. GVCF written to ../../${FILENAME}.g.vcf" >> ${LOG_FILE} 

}


mkdir -p ${LOG_PATH} ${OUTPATH} ${BAM_OUT}
export -f run_haplotypeCaller

find ${INPATH} -maxdepth 1 -type f -name '*.bam' | parallel -j 5 --progress --verbose run_haplotypeCaller {}

