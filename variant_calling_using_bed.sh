#!/bin/bash
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --job-name=bcftools_window_call
#SBATCH --output=logs/bcftools_window_%A_%a.out
#SBATCH --error=logs/bcftools_window_%A_%a.err
#SBATCH --array=1-21
#SBATCH --cpus-per-task=64
#SBATCH --mem=450G
#SBATCH --time=96:00:00
#SBATCH --tmp=300G

# Set paths
WORKING_DIR="/scratch.global/qiuxx221/ebi_mapping"
BED_FILE="${WORKING_DIR}/bed_windows/barley_3x_windows.bed"
BAM_LIST="${WORKING_DIR}/BAM_dedup/done/276_bam_list.txt"
REF="/home/steffenb/shared/ris/barley_genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa"
OUT_DIR="${WORKING_DIR}/raw-VCF"
mkdir -p "$OUT_DIR" logs

# Extract the region line (1-based for SLURM)
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BED_FILE")

# Parse region
CHR=$(echo "$LINE" | cut -f1)
START=$(echo "$LINE" | cut -f2)
END=$(echo "$LINE" | cut -f3)

# Create region string
REGION="${CHR}:${START}-${END}"
OUT_VCF="${OUT_DIR}/region_${CHR}_${START}_${END}.vcf"

# Run bcftools on this region
bcftools mpileup \
  -Ou \
  --threads 64 \
  -f "$REF" \
  --bam-list "$BAM_LIST" \
  -r "$REGION" \
  -q 5 \
  -I \
  -a FMT/AD,FMT/DP | \
bcftools call \
  --threads 64 \
  -mv \
  -f GQ \
  -Ov \
  > "$OUT_VCF"
