#!/bin/bash
#
# Script for performing QC of raw sequencing data

# Specify paths
FASTQ_PATH=/node1nfs/archive/data/trifolium/toronto_gwsd/128.120.88.251/X202SC20030702-Z01-F001/raw_data/

# Outpath depends on whether QC is being performed on 
# raw FASTQ files or those post-adapter trimming
OUTPATH=../../analysis/sequencingQC/post_trimming_reports/
#OUTPATH=../../analysis/sequencingQC/raw_fastq_reports/

# Run FASTQC in parallel on all raw sequence files
# find '*_1.fq.gz' if raw sequences, '*_trimmed.fq.gz' if post trimming.
find ${FASTQ_PATH} -name '*_trimmed.fq.gz' | parallel -j 30 --verbose \
 fastqc {} --outdir ${OUTPATH}
