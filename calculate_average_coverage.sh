#!/usr/bin/env bash
#SBATCH --job-name=bam_depth
#SBATCH --output=logs/depth_%A_%a.out
#SBATCH --error=logs/depth_%A_%a.err
#SBATCH --array=0-276    
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module load samtools
# Set BAM directory
BAM_DIR="/scratch.global/qiuxx221/ebi_mapping/BAM_dedup"
OUT_DIR="${BAM_DIR}/depth_results"
mkdir -p "$OUT_DIR"

# List BAMs into an array
BAM_LIST=($(ls "$BAM_DIR"/*.bam))

# Select BAM based on SLURM_ARRAY_TASK_ID
BAM="${BAM_LIST[$SLURM_ARRAY_TASK_ID]}"
BAM_BASENAME=$(basename "$BAM" .bam)
OUTFILE="${OUT_DIR}/${BAM_BASENAME}_avg_depth.txt"

# Index BAM if not indexed
if [ ! -f "${BAM}.csi" ] && [ ! -f "${BAM}.bai" ]; then
    samtools index -c "$BAM"
fi

# Calculate average depth
avg=$(samtools depth "$BAM" | awk '{sum+=$3} END { if (NR>0) printf "%.2f", sum/NR; else print "NA" }')

# Write output
printf "BAM_File\tAverage_Depth\n%s\t%s\n" "$(basename "$BAM")" "$avg" > "$OUTFILE"

echo "Done. Result saved to $OUTFILE"
