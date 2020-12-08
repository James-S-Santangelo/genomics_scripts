#!/bin/bash
#
# Script to QC raw BAM files using Qualimap

# Required to avoid X11 error
unset DISPLAY 

export LOG_FILE=./log_files/qualimap/rawBamQC_Qualimap.log
export BAM_PATH=../../data/raw/bam/bwa_mapping
export OUTPATH=../../analysis/mappingQC/raw_bam_qc

{
~/Downloads/qualimap_v2.2.1/qualimap multi-bamqc \
    --paint-chromosome-limits \
    --run-bamqc \
    --data ../../resources/bamPaths_forQualimap.txt \
    -outdir ${OUTPATH} \
    --java-mem-size=300G &&

# Need to manually  move bamQC runs from qualimap to multi-mabQC output directory
mkdir -p ${OUTPATH}/bamQC_perSample
mv ${BAM_PATH}/*_stats ${OUTPATH}/bamQC_perSample
} > ${LOG_FILE}

