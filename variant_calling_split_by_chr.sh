#!/bin/bash
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --job-name=bcftools_call
#SBATCH --output=logs/bcftools_call_%A_%a.out
#SBATCH --error=logs/bcftools_call_%A_%a.err
#SBATCH --array=1-7
#SBATCH --cpus-per-task=64
#SBATCH --mem=450G
#SBATCH --time=96:00:00
#SBATCH --tmp=300G

# Load modules (modify for your system)
# module load bcftools

# Define working directories and inputs
WORKING_DIR="/scratch.global/qiuxx221/ebi_mapping"
REF="/home/steffenb/shared/ris/barley_genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa"
BAM_LIST="/scratch.global/qiuxx221/ebi_mapping/BAM_dedup/done/276_bam_list.txt"
OUT_DIR="${WORKING_DIR}/raw-VCF"

mkdir ${WORKING_DIR}/raw-VCF

# Chromosome list (edit as needed)
CHR_LIST=(1H 2H 3H 4H 5H 6H 7H)
CHROM=${CHR_LIST[$SLURM_ARRAY_TASK_ID-1]}

# Output file
OUT_VCF="${OUT_DIR}/${CHROM}-subpop.vcf"

# Create output and log directories if not exist
mkdir -p "$OUT_DIR"
mkdir -p logs

# Run bcftools pipeline
bcftools mpileup \
   -Ou \
   --threads 64 \
   -f "$REF" \
   --bam-list "$BAM_LIST" \
   -q 5 \
   -r "$CHROM" \
   -I \
   -a FMT/AD,FMT/DP | \
bcftools call \
   --threads 64 \
   -mv \
   -f GQ \
   -Ov \
   > "$OUT_VCF"
