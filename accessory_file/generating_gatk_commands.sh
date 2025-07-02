#!/bin/bash

REF="/scratch.global/qiuxx221/ebi_mapping/barley_genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa"
BAM_LIST="/scratch.global/qiuxx221/ebi_mapping/RG/filter_mapq20/281_bam_list_mapq20.txt"
INTERVALS="/scratch.global/qiuxx221/ebi_mapping/barley_genome/genome_2Mb_intervals.list"
OUTDIR="/home/steffenb/public/qiuxx221/gatk_joint_output_mapq20_parallel"
mkdir -p "$OUTDIR"

# === OUTPUT FILE ===
CMD_FILE="gatk_joint_commands_para.txt"
> "$CMD_FILE"

# === Build a single BAM argument string with multiple -I flags ===
INPUT_ARGS=""
while read -r BAM; do
  INPUT_ARGS+=" -I $BAM"
done < "$BAM_LIST"

# === Loop over intervals and build commands ===
while read -r INTERVAL; do
  REGION=$(echo "$INTERVAL" | tr ':-' '_')  # e.g. 3H_42000001_43000000
  OUTVCF="${OUTDIR}/${REGION}.mapq20.vcf.gz"

  echo "gatk --java-options \"-Xmx40g -Djava.io.tmpdir=/scratch.global/\$SLURM_JOB_ID\" HaplotypeCaller -R $REF $INPUT_ARGS -L $INTERVAL -O $OUTVCF ; tabix -p vcf $OU
TVCF " >> "$CMD_FILE"
done < "$INTERVALS"
