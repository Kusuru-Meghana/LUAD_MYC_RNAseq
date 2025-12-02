#!/bin/bash

# ------------------------------
# LAYER 2 RNA-seq - Step 1
# Download SRA files
# Convert SRA → FASTQ (single-end)
# ------------------------------

# Move to raw_fastq folder
cd /mnt/c/Users/megha/projects/MYC_Project/LAYER2_RNAseq/raw_fastq

# Download SRA files
prefetch SRR8309421
prefetch SRR8309422
prefetch SRR8309423
prefetch SRR8309424

# Convert SRA → FASTQ (single-end)
fasterq-dump SRR8309421/SRR8309421.sra --split-files -O . -e 8
fasterq-dump SRR8309422/SRR8309422.sra --split-files -O . -e 8
fasterq-dump SRR8309423/SRR8309423.sra --split-files -O . -e 8
fasterq-dump SRR8309424/SRR8309424.sra --split-files -O . -e 8

# Compress FASTQ files
gzip SRR8309421.fastq
gzip SRR8309422.fastq
gzip SRR8309423.fastq
gzip SRR8309424.fastq
