# Differential Expression Results (DESeq2 Output)

This folder contains all differential expression results generated during **Layer 2 RNA-seq analysis** of A549 MYC knockdown (siMYC) and control samples.

Differential expression analysis was performed using **DESeq2** in R.

---

## 1. Files in this folder

### **DESeq2_siMYC_vs_CTRL_results.xlsx**
- Complete DESeq2 results table.
- Includes:
  - log2 fold change  
  - p-values  
  - adjusted p-values (FDR)  
  - baseMean (expression strength)  
  - gene identifiers  

### **differential_expression.R**
- The R script used to run the full differential expression pipeline.
- Contains:
  - Import of featureCounts matrix  
  - DESeq2 normalization  
  - PCA  
  - MA plot  
  - Volcano plot  
  - Heatmap of top variable genes  
  - Export of results for enrichment analysis  

Direct link to the script:
- [`differential_expression.R`](./differential_expression)

---

## 2. Summary of Key Findings

### **Downregulated in siMYC**
Genes dependent on MYC activity showed decreased expression:
- RNA splicing factors  
- Ribosome biogenesis genes  
- DNA replication machinery  
- Chromatin assembly components  
- Translational regulators  

These represent classic MYC-driven growth programs.

### **Upregulated in siMYC**
A strong interferon/immune signature was activated:
- ISG15  
- IFIT1, IFIT2, IFIT3  
- IRF1  
- STAT1  
- OAS1, OAS2, OASL  
- MX1  
- BST2  
- CMPK2  

These genes become derepressed when MYC levels are reduced.

---

## 3. Figures Generated from This Analysis

All plots generated from DESeq2 are stored in the `plots/` folder:

- PCA plot  
- MA plot  
- Volcano plot  
- Heatmap (Top 50 DEGs)  

These visualizations provide quality assessment and confirm distinct clustering between CTRL and siMYC groups.

---

## 4. Purpose

The DESeq2 results in this folder are used for:

- Visualization of sample separation  
- Identification of MYC-regulated genes  
- Functional enrichment analysis (GO, KEGG, Reactome)  
- Integration with MYC ChIP-seq data (Layer 1)  
- Construction of Layer 2 biological conclusions  

---

## 5. Notes

- DESeq2 normalization automatically adjusts for sequencing depth and dispersion.  
- Fold changes reflect transcriptional reprogramming following MYC knockdown.  
- This repository only includes processed DESeq2 results; raw FASTQ and BAM files are excluded due to size limitations.

