#!/bin/bash
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --job-name=bcfsubset
#SBATCH --output=bcfconcat.out
#SBATCH --error=bcfconcat.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=48:00:00
#SBATCH --tmp=300G

# this file is used to join 14 vcf files together as a whole
module load bcftools 
bcftools concat -f vcf_merge_list.txt -Oz -o bcftools_repadapt_mapq20_concatenated.vcf.gz

#bcftools index concatenated.vcf.gz
bcftools index --csi bcftools_repadapt_mapq20_concatenated.vcf.gz
