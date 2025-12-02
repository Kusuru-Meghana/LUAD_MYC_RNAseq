#!/bin/bash

# --------------------------------------------------
# LAYER 2 RNA-seq - Step 2
# Quality control (FastQC + MultiQC)
# --------------------------------------------------

# Move to qc folder
cd /mnt/c/Users/megha/projects/MYC_Project/LAYER2_RNAseq/qc

# Run FastQC on all FASTQ files
fastqc ../raw_fastq/*.fastq.gz -o .

# Run MultiQC to aggregate reports
multiqc . -o .
