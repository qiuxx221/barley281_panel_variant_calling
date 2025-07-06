#!/bin/bash
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --job-name=merge_indels
#SBATCH --output=merge_indels%A_%a.out
#SBATCH --error=merge_indels%A_%a.err
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

bcftools concat \
    1H.concat_indels.vcf.recode.vcf \
    2H.concat_indels.vcf.recode.vcf \
    3H.concat_indels.vcf.recode.vcf \
    4H.concat_indels.vcf.recode.vcf \
    5H.concat_indels.vcf.recode.vcf \
    6H.concat_indels.vcf.recode.vcf \
    7H.concat_indels.vcf.recode.vcf \
    -Ov -o all_chromosomes.gatk_mapq20_raw_indels.vcf

bgzip all_chromosomes.gatk_mapq20_raw_indels.vcf
bcftools index -c all_chromosomes.gatk_mapq20_raw_indels.vcf.gz
