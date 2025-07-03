#!/bin/bash
#SBATCH --job-name=intersect_callable
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --output=intersect_callable_%A_%a.out
#SBATCH --error=intersect_callable_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --array=1-14
#SBATCH --mail-user=qiuxx221@umn.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load bcftools/1.21

# Input directory with VCFs
vcf_dir="/scratch.global/qiuxx221/ebi_mapping/mapq20-VCF/split_by_14/new_filter"
# BED file to intersect with
bed_file="/home/steffenb/shared/ris/barley281_variant_call/callable/callable_bed_convert"
# Output directory
outdir="${vcf_dir}/intersected"
mkdir -p "$outdir"

# File list with VCF filenames (one per line)
file_list="${vcf_dir}/14_tag_gz_files.txt"

# Get filename for this array index (1-based)
vcf_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$file_list")
base_name=$(basename "$vcf_file" .tags.vcf.gz)

# Run intersection
bcftools view -R "$bed_file" "${vcf_dir}/${vcf_file}" -Oz -o "${outdir}/${base_name}.intersect.vcf.gz"
bcftools index "${outdir}/${base_name}.intersect.vcf.gz"
