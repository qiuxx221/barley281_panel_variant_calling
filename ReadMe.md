The Variant Calling for the 281 wild barley collection was performed using Repadapt workflow as the foundation (https://github.com/JimWhiting91/RepAdapt/tree/main/snp_calling_pipeline). 

# Reads trimming and Mapping

raw reads were first trimmed using fastp program using and aligned to the morex v3 reference genome via bwa mem (https://github.com/SteffensonLab/Barley_IPK_variant_calling/blob/main/01_snp_calling-bash/03_map.sh) ). Picard tool (version 2.25.6) was used to mark duplicates. Duplicated reads were removed from the bam file (https://github.com/SteffensonLab/Barley_IPK_variant_calling/blob/main/01_snp_calling-bash/04_remove_dup.sh) and the final bam file was filtered using mapping quality of 20 (https://github.com/qiuxx221/barley281_panel/blob/main/mapq20_filters_bam.sh)

# Variant calling
Two methods were used for the variant calling: bcftools mpileup and GATK

## Bcftools mpileup variant calling

Due to the large barley chromosome size, each chromosome were split in half to speed up the variant calling process, this result in 14 calling windows (https://github.com/qiuxx221/barley281_panel/blob/main/accessory_file/all_7_chr_2x.region).

281 genotype variant call was done using the panel joint calling (https://github.com/qiuxx221/barley281_panel/blob/main/bcftools_mpileup_variant_calling.sh). 

To ensure the high quality set of SNPs were obtained, Variant was further hard filtered using the miniumn read depth of 10, and quality of 30 (https://github.com/qiuxx221/barley281_panel/blob/main/snp_filtering.sh). After the quality filtering by site and by genotype, sites with no variants were removed. Allele frequency, allele count, allele number were recalculated using bcftools +fill-tags function. The final vcf file was created by concatenating the 14 vcfs together.

## GATK variant calling

Limited by the speed of gatk haplotype caller, barley genome was divided into 2 million base pair window (https://github.com/qiuxx221/barley281_panel/blob/main/accessory_file/genome_2Mb_intervals.list) for the joint calling using bam files filtered by mapping quality of 20 . This result in a total of 2102 windows

Script generating_gatk_commands.sh was used to generate the variant calling command https://github.com/qiuxx221/barley281_panel/blob/main/accessory_file/generating_gatk_commands.sh. Those commands were further split into 24 commands to run parallel. 

One tricky part of the gatk calling is that there are issuing when the chromomsome positons are after 536,000,000 bp. So for the final variant calls, we first run standard gatk processing using commands at https://github.com/qiuxx221/barley281_panel/blob/main/accessory_file/gatk_commands.zip, then we run https://github.com/qiuxx221/barley281_panel/blob/main/accessory_file/csi_window_rerun_cmd.zip to generate those out of bound window variant call. 

Finally, results were combined and cleaned to ensure all windows were coveraged properly and merged together. 

To ensure the high quality variant calls, standard gatk filtering paramater was used. 

# Variant annotation

To annotate the SNPs, variant effector predictor was used using commands at (https://github.com/qiuxx221/barley281_panel/blob/main/vep_SNP_annotation.sh)

The conda environment was exported as a yaml file and stored at https://github.com/qiuxx221/barley281_panel/blob/main/accessory_file/vep_environment.yaml

sorted gff files for this annotation was stored at https://github.com/qiuxx221/barley281_panel/blob/main/accessory_file/Hordeum_vulgare.MorexV3.sorted.gff3.gz alone with its csi index 





