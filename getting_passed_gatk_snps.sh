#!/bin/bash
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --job-name=filtered_snp_merge
#SBATCH --output=filtered_snp_merge%A_%a.out
#SBATCH --error=filtered_snp_merge%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=qiuxx221@umn.edu
#SBATCH --tmp=300Gb

# Load environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate apg_migrate
module load bcftools/1.21
module load gatk


bgzip -@8 1H.concat_filtered-snps.vcf
bcftools index -c 1H.concat_filtered-snps.vcf
# Input VCF list file
vcf_list="filtered_snps_list.txt"

# Concat all VCFs listed
bcftools concat -f "$vcf_list" -Ov -o all_chromosomes.gatk_mapq20_filtered_snp.vcf

# Compress and index
bgzip -@8 -f  all_chromosomes.gatk_mapq20_filtered_snp.vcf
bcftools index -c  all_chromosomes.gatk_mapq20_filtered_snp.vcf.gz

bcftools view -f PASS -Oz -o all_chromosomes.gatk_mapq20_filtered_snp.pass.vcf.gz all_chromosomes.gatk_mapq20_filtered_snp.vcf.gz
