#!/bin/bash
#
# Script to perform Qualimap QC on BAMs foolowing mapping or trimmed reads using either BWA or NGM

unset DISPLAY # Required to avoid X11 error

~/Downloads/qualimap_v2.2.1/qualimap multi-bamqc -c \
    --run-bamqc \
    --data ../../resources/qualimap_bwa-ngm_compare.config \
    -outdir ../../analysis/mappingQC/compare_mapping \
    --java-mem-size=100G
