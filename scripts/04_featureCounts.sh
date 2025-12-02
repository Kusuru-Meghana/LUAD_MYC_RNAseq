#!/bin/bash

# ===========================================
# LAYER 2 RNA-seq : FeatureCounts Script
# Count reads per gene for 4 samples
# Single-end RNA-seq
# ===========================================

# ---- Paths ----
BAM_DIR="../aligned_bam"
GTF="../genome_index/gencode.v45.annotation.gtf"
OUT_DIR="../counts"

# Create output folder if missing
mkdir -p $OUT_DIR

# ---- Run featureCounts ----
featureCounts \
    -T 8 \
    -s 0 \
    -a $GTF \
    -o $OUT_DIR/gene_counts.txt \
    $BAM_DIR/SRR8309421.sorted.bam \
    $BAM_DIR/SRR8309422.sorted.bam \
    $BAM_DIR/SRR8309423.sorted.bam \
    $BAM_DIR/SRR8309424.sorted.bam

echo " FeatureCounts Finished!"
echo "Output files created in: $OUT_DIR/"
