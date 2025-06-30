#!/bin/bash
#SBATCH --job-name=annotate_SNPs_all
#SBATCH --output=annotate_SNPs_all
#SBATCH --error=annotate_SNPs_all
#SBATCH --time=12:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

source $(conda info --base)/etc/profile.d/conda.sh
conda activate ensembl_env_reinstall

/users/5/qiuxx221/.conda/envs/ensembl_env_reinstall/bin/vep \
    -i /home/steffenb/shared/ris/barley281_variant_call/repadapt_mapq20_vcf/1final_info_genotype_filtered.vcf.gz \
    --gff /scratch.global/qiuxx221/ebi_mapping/barley_genome/Hordeum_vulgare.MorexV3.sorted.gff3.gz \
    --fasta /scratch.global/qiuxx221/ebi_mapping/barley_genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa \
    --format vcf \
    --fork 128 \
    --force_overwrite \
    --output_file repadapt_final_info_genotype_filtered.vcf.annotation.txt


# the annotation has now been moved to /home/steffenb/shared/ris/barley_genome/annotation_vep to host permanently 
