# Functional Enrichment Analysis (GO, KEGG, Reactome)

This folder contains all pathway and functional enrichment results generated during **Layer 2 RNA-seq analysis** of A549 MYC knockdown (siMYC) vs control samples.

Enrichment analysis was performed using the R package **clusterProfiler**, along with **org.Hs.eg.db**, **ReactomePA**, and KEGG pathway databases.

---

## 1. Files in this folder

### **GO_Enrichment.xlsx**
- Gene Ontology enrichment results for:
  - Biological Process (BP)
  - Molecular Function (MF)
  - Cellular Component (CC)
- Column details:
  - Description  
  - GeneRatio  
  - p.adjust  
  - Count  
  - Enriched gene symbols  

### **KEGG_Enrichment.xlsx**
- KEGG pathway enrichment results.
- Includes:
  - Significantly upregulated pathways  
  - Pathways downregulated following MYC knockdown  

### **Reactome_Enrichment.xlsx**
- Reactome pathway enrichment results.
- Contains:
  - Pathway names  
  - Enrichment score  
  - FDR  
  - Gene lists  

### **enrichment_analysis.R**
- R script used to generate all enrichment results.
- Includes:
  - ID conversion (ENSG → SYMBOL)  
  - GO (BP/MF/CC) enrichment  
  - KEGG pathway analysis  
  - Reactome pathway analysis  
  - Export of tables for visualization  

Direct link:
- [`enrichment_analysis.R`](../scripts/enrichment_analysis)

---

## 2. Summary of Enrichment Results

### **Downregulated pathways (MYC-activated programs)**  
These pathways decreased after MYC knockdown and represent classical MYC-driven biological functions:

- **RNA splicing & mRNA processing**  
- **Ribosome biogenesis**  
- **Nucleosome assembly & chromatin organization**  
- **DNA replication & cell cycle progression**  
- **RNA polymerase I transcription**

These findings reflect MYC’s role as a global amplifier of transcription and growth-related processes in cancer.

---

### **Upregulated pathways (MYC-repressed programs)**  
An interferon/innate immune signature becomes strongly enriched:

- **Interferon signaling**  
- **Innate immune activation**  
- **Viral defense pathways**  
- **Antigen presentation**  

Genes driving this signal include:

- ISG15  
- IFIT1/2/3  
- OAS1/2/OASL  
- IRF1  
- STAT1  
- MX1  
- BST2  

These pathways are typically **suppressed by MYC**, becoming derepressed when MYC levels are lowered.

---

## 3. Interpretation

The enrichment results highlight the **dual role of MYC**:

### **1. MYC activates**
- Transcription  
- RNA processing  
- Ribosome formation  
- Cell cycle  
- Chromatin assembly  

Loss of MYC leads to widespread downregulation of these programs.

### **2. MYC represses**
- Interferon-stimulated genes  
- Immune surveillance pathways  
- Innate antiviral response  

These pathways become activated in the siMYC condition.

This represents a classical MYC-loss expression phenotype observed in LUAD and other MYC-driven cancers.

---

## 4. Purpose of These Results

The enrichment results support:

- Biological interpretation of DEGs  
- Construction of the Layer 2 MYC regulatory map  
- Pathway-level integration with MYC ChIP-seq (Layer 1)
- Visualization in the upcoming multi-omics (Layer 1 + Layer 2) repository  

---

## 5. Notes

- All enrichment analyses were performed using gene symbols.  
- Adjusted p-values (FDR) were used to determine significance.  
- Raw FASTQ/BAM files are excluded from this repository due to size limitations.

