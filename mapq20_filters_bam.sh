#!/bin/bash
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --job-name=mapq_filter
#SBATCH --output=logs/mapq_filter_%A_%a.out
#SBATCH --error=logs/mapq_filter_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=1-281 # Replace <N> with number of BAM files
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qiuxx221

# Load samtools if needed
module load samtools

# Input/output directories
IN_DIR="/scratch.global/qiuxx221/ebi_mapping/RG"
OUT_DIR="/scratch.global/qiuxx221/ebi_mapping/RG/filter_mapq20"
mkdir -p "$OUT_DIR"

# Get list of BAM files
BAM_LIST=($(ls "$IN_DIR"/*.bam))
BAM="${BAM_LIST[$SLURM_ARRAY_TASK_ID-1]}"
BASENAME=$(basename "$BAM" .bam)
OUT_BAM="${OUT_DIR}/${BASENAME}_MAPQ20.bam"

echo "[$SLURM_ARRAY_TASK_ID] Processing $BASENAME.bam..."

# Filter and index
samtools view -h -q 20 -b "$BAM" > "$OUT_BAM"
samtools index -c "$OUT_BAM"

echo "[$SLURM_ARRAY_TASK_ID] Done with $OUT_BAM"

