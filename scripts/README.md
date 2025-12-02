# RNA-seq Processing Scripts

This folder contains all shell scripts used for the **Layer 2 RNA-seq workflow** of A549 MYC knockdown (siMYC) analysis.

Each script corresponds to one step of the RNA-seq pipeline, from FASTQ download to read counting.  
The scripts are designed to run sequentially.

---

## 1. `01_download_fastq.sh`

**Purpose:**  
Downloads raw FASTQ files from SRA and converts `.sra` → `.fastq.gz`.

**Includes:**  
- `prefetch` for downloading SRA files  
- `fasterq-dump` for extracting FASTQ  
- gzip compression  

**Script:**  
[`01_download_fastq.sh`](01_download_fastq.sh)

---

## 2. `02_qc.sh`

**Purpose:**  
Runs quality control on all FASTQ files.

**Includes:**  
- FastQC for per-sample QC  
- MultiQC for aggregated summary  

**Script:**  
[`02_qc.sh`](02_qc.sh)

---

## 3. `03_align_hisat2.sh`

**Purpose:**  
Aligns reads to the human genome (GENCODE hg38) using HISAT2.

**Includes:**  
- Single-end HISAT2 alignment  
- SAM → BAM conversion  
- Sorted BAM creation  
- BAM indexing  
- Alignment logs  

**Script:**  
[`03_align_hisat2.sh`](03_align_hisat2.sh)

---

## 4. `04_featureCounts.sh`

**Purpose:**  
Generates gene-level read counts using featureCounts.

**Includes:**  
- GENCODE v45 annotation  
- Single-end mode  
- Outputs `gene_counts.txt` for DESeq2 analysis  

**Script:**  
[`04_featureCounts.sh`](04_featureCounts.sh)

---

## Notes

- Scripts assume project structure:
