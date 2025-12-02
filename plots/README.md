# RNA-seq Visualization Plots

This folder contains all visual outputs generated during the **Layer 2 RNA-seq analysis** of A549 MYC knockdown (siMYC) and control samples.  
Plots include quality assessment, differential expression visualizations, sample clustering, and functional enrichment.

---

# 1. Dispersion Plot (DESeq2)

**File:** `Dispersion_plot.png`

This plot shows the dispersion estimates for each gene after DESeq2 modeling.

- Black: gene-wise dispersion estimates  
- Red: fitted dispersion trend  
- Blue: final dispersion values after shrinkage  

**Interpretation:**  
A clean, smooth red trend confirms correct variance estimation, validating the dataset for reliable differential expression analysis.

---

# 2. MA Plot

**File:** `MA_plot.png`

Shows mean expression vs log2 fold change.

- Grey: non-significant genes  
- Blue: significantly differentially expressed genes  

**Interpretation:**  
MYC knockdown causes broad downregulation of highly expressed genes, consistent with MYC acting as a transcriptional amplifier.

---

# 3. PCA Plot

**File:** `PCA_plot.png`

Principal Component Analysis of normalized RNA-seq counts.

- PC1 = 91% of variance  
- CTRL and siMYC samples form two distinct clusters  

**Interpretation:**  
The separation indicates MYC knockdown is the primary source of variation, showing a strong biological effect with no major batch effects.

---

# 4. Heatmap (Top 50 DEGs)

**File:** `Heatmap_top50_DE_genes.png`

Displays the top 50 most variable differentially expressed genes.

- Blue = downregulated  
- Red = upregulated  
- Hierarchical clustering shows perfect separation of CTRL vs siMYC  

**Interpretation:**  
Consistent expression patterns confirm the robustness of biological replicates.

---

# 5. GO Biological Process (BP) Enrichment — Dot Plot

**File:** `GO_BP_dotplot.png`

Displays enriched **Biological Process** GO terms for all DEGs.

Most enriched categories include:

- RNA splicing (multiple terms)  
- Protein–DNA complex assembly  
- Nucleosome assembly  
- Regulation of mRNA processing  

**Interpretation:**  
MYC knockdown strongly impacts transcriptional and RNA processing pathways — core MYC functions in cancer cells.

---

# 6. KEGG Pathway Enrichment — Dot Plot

**File:** `KEGG_dotplot.png`

Top enriched KEGG pathways include:

- Cell cycle (down)  
- p53 signaling (up)  
- Cellular senescence  
- Chromatin remodeling  
- Interferon-related viral pathways (up)  

**Interpretation:**  
Loss of MYC activity induces stress and immune-response pathways while suppressing cell proliferation pathways.

---

# 7. Reactome Pathway Enrichment — Dot Plot

**File:** `Reactome_dotplot.png`

Highly enriched Reactome pathways include:

- DNA replication  
- Senescence-associated secretory phenotype (SASP)  
- rRNA expression and ribosome biogenesis  
- RNA Polymerase I promoter opening  
- Chromosome condensation  

**Interpretation:**  
Reactome confirms MYC’s central role in cell cycle progression and rRNA production.  
MYC knockdown leads to strong senescence-like and stress responses.

---

# Summary

These plots demonstrate:

- Excellent sequencing quality  
- Strong transcriptional reprogramming upon MYC knockdown  
- Distinct clustering between CTRL and siMYC  
- Downregulation of proliferative pathways  
- Upregulation of immune and stress pathways  

Together, they provide visual confirmation of MYC’s biological roles in LUAD cells and support all downstream interpretations in the Layer 2 RNA-seq analysis.

