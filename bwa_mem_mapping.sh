#!/bin/bash

# This assumes SAMPLE is passed in as an environment variable
THREADS=8
INPUT_DIR="/scratch.global/qiuxx221/ebi_mapping/reads"
WORKING_DIR="/scratch.global/qiuxx221/ebi_mapping/bwa_output"
REF="/home/steffenb/shared/ris/barley_genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa"

# Input files
R1="${INPUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
R2="${INPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
OUT_BAM="${WORKING_DIR}/${SAMPLE}.bam"
SORTED_BAM="${WORKING_DIR}/${SAMPLE}_sorted.bam"

# Ensure output dir exists
mkdir -p "${WORKING_DIR}"

# Run pipeline
bwa mem -t $THREADS "$REF" "$R1" "$R2" | \
samtools view -Sb -q 10 -@ $THREADS > "$OUT_BAM" && \
samtools sort -@ $THREADS -o "$SORTED_BAM" "$OUT_BAM" && \
samtools index --csi "$SORTED_BAM" && \
rm -f "$OUT_BAM"
