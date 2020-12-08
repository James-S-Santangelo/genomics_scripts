#!/bin/bash
#
# Script to compare read mapping to the reference using both BWA and NextGenMap

INPATH=/node1nfs/archive/data/trifolium/toronto_gwsd/128.120.88.251/X202SC20030702-Z01-F001/raw_data/
samples=(s_100_13 s_116_1 s_23_12 s_37_18 s_40_8 s_43_15 s_54_9 s_7_20 s_83_10 s_97_3)

export BWA_OUTPATH=/scratch/research/projects/trifolium/GWSD_Genome-Wide-Delection-Demography/data/raw/bam/compare_mapping/bwa
export NGM_OUTPATH=/scratch/research/projects/trifolium/GWSD_Genome-Wide-Delection-Demography/data/raw/bam/compare_mapping/ngm

export PATH_TO_REFERENCE=/scratch/research/references/trifolium/repens/GCA_005869975.1_AgR_To_v5/GCA_005869975.1_AgR_To_v5_genomic.fna.gz

fullPathsToSamples=()
for ((i = 0; i < ${#samples[@]}; i++))
do
    fullPathsToSamples[$i]=${INPATH}${samples[$i]} 
    #echo "${amples[$i]}"
done

function run_bwa {
    
    # Inputs and outputs
    local SAMPLE=$(basename ${1})
    local READ_1=$(find ${1} -type f -name '*_1_trimmed.fq.gz' -printf "%f\n")
    local READ_2=$(find ${1} -type f -name '*_2_trimmed.fq.gz' -printf "%f\n")
    local PATH_TO_READ_1="${1}/trimmed/${READ_1}"
    local PATH_TO_READ_2="${1}/trimmed/${READ_2}"
    local OUTFILE=${BWA_OUTPATH}/${SAMPLE} 
    
    echo "Processing reads for sample ${SAMPLE}"
    sleep 1

    echo "Aligning reads for ${READ_1} and ${REAND_2} to reference genome using bwa. Storing mapping as ${OUTFILE}.sam"
    sleep 1
    bwa mem ${PATH_TO_REFERENCE} ${PATH_TO_READ_1} ${PATH_TO_READ_2} > "${OUTFILE}.sam"  

    echo "Converting ${OUTFILE}.sam to bamfile ${OUTFILE}.bam" 
    sleep 1
    samtools view -b "${OUTFILE}.sam" > "${OUTFILE}.bam" &&

    echo "Sorting ${OUTFILE}.bam"
    sleep 1
    samtools sort "${OUTFILE}.bam" > "${OUTFILE}_sorted.bam" &&

    echo "Indexing "${OUTFILE}_sorted.bam"
    sleep 1
    samtools index ${OUTFILE}_sorted.bam" &&

    echo "Done mapping reads for ${SAMPLE}"
    sleep 1

    echo "Removing temporary files"
    rm "${OUTFILE}.bam"
    rm "${OUTFILE}.sam"
}


function run_ngm {
    
    # Inputs and outputs
    local SAMPLE=$(basename ${1})
    local READ_1=$(find ${1} -type f -name '*_1_trimmed.fq.gz' -printf "%f\n")
    local READ_2=$(find ${1} -type f -name '*_2_trimmed.fq.gz' -printf "%f\n")
    local PATH_TO_READ_1="${1}/trimmed/${READ_1}"
    local PATH_TO_READ_2="${1}/trimmed/${READ_2}"
    local INTERLEAVED_OUT="${1}/trimmed/${SAMPLE}_interleaved.fq"
    local INTERLEAVED_BGZIP="${INTERLEAVED_OUT}.gz"
    local OUTFILE=${NGM_OUTPATH}/${SAMPLE} 
    
    echo "Processing reads for sample ${SAMPLE}"
    sleep 1

    # NextGenMap requires interleved reads if forward and reverse reads are not 
    # in the exact same order, which is likely if we lost reads post-trimming.
    echo "Interleaving reads for ${SAMPLE}. Suggested for NGM since it doesn't check read names"
    ngm-utils interleave -1 ${PATH_TO_READ_1} -2 ${PATH_TO_READ_2} -o ${INTERLEAVED_OUT} &&

    bgzip --threads 2 --force ${INTERLEAVED_OUT} > ${INTERLEAVED_BGZIP} &&
    
    echo "Aligning reads for ${READ_1} and ${READ_2} to reference genome using ngm. Storing mapping as ${OUTFILE}.sam"
    sleep 1
    
    ngm --ref ${PATH_TO_REFERENCE} \
        --qry ${INTERLEAVED_BGZIP} \
        --paired \
        --output "${OUTFILE}.sam" && 

    echo "Converting ${OUTFILE}.sam to bamfile ${OUTFILE}.bam" 
    sleep 1
    samtools view -b "${OUTFILE}.sam" > "${OUTFILE}.bam" &&

    echo "Sorting ${OUTFILE}.bam"
    sleep 1
    samtools sort "${OUTFILE}.bam" > "${OUTFILE}_sorted.bam" &&

    echo "Indexing "${OUTFILE}_sorted.bam"
    sleep 1
    samtools index ${OUTFILE}_sorted.bam" &&

    echo "Done mapping reads for ${SAMPLE}"
    sleep 1

    echo "Removing temporary files"
    rm "${INTERLEAVE_OUT}"
    rm "${INTERLEAVED_BGZIP}"
    rm "${OUTFILE}.bam"
    rm "${OUTFILE}.sam"
}

export -f run_bwa
export -f run_ngm

# Run bwa and ngm in parallel 
parallel -j 10 --progress --verbose run_bwa ::: ${fullPathsToSamples[@]} & 
parallel -j 10 --progress --verbose run_ngm ::: ${fullPathsToSamples[@]} 
