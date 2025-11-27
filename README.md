# A549 MYC Knockdown (siMYC) RNA-seq Analysis — Layer 2 of LUAD MYC Multi-Omics

This repository contains a complete RNA-seq workflow for A549 lung adenocarcinoma cells following MYC knockdown using siRNA (siMYC).
The objective is to characterize transcriptomic changes driven by loss of MYC, a major oncogenic transcription factor in LUAD.

Although this analysis is part of a larger multi-omics framework, this repository is fully self-contained and functions as a stand-alone RNA-seq project.
It represents Layer 2 of the LUAD MYC regulatory map:

- Layer 1: MYC ChIP-seq (DNA binding sites)

- Layer 2: RNA-seq (this repo — MYC-dependent gene expression)

- Layer 3: ATAC-seq (chromatin accessibility, coming next)

A separate combined repository will later hold cross-layer figures (Layer1+Layer2).

## 1. Introduction

MYC is a master regulator of cell growth, metabolism, ribosome biogenesis, and transcriptional amplification in cancer.
To understand how MYC controls gene expression in LUAD, the MYC gene was silenced using siRNA, generating:

- CTRL: A549 cells with normal MYC

- siMYC: A549 cells with reduced MYC expression

RNA-seq enables the identification of:

- Genes that require MYC for expression (downregulated in siMYC)

- Genes that MYC normally represses (upregulated in siMYC)

- Pathways impacted by MYC loss

This provides a quantitative map of MYC-dependent transcriptional programs.


## 2. Dataset Overview
| Sample | Condition | SRA Accession |
| ------ | --------- | ------------- |
| CTRL1  | Control   | SRR8309421    |
| CTRL2  | Control   | SRR8309422    |
| siMYC1 | Knockdown | SRR8309423    |
| siMYC2 | Knockdown | SRR8309424    |
All processing in this project uses raw FASTQ files downloaded from SRA.

## 3. Workflow Summary
