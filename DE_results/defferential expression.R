setwd("C:/Users/megha/projects/MYC_Project/LAYER2_RNAseq")
getwd()   # just to check

#Install & load DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")



install.packages(c("tidyverse", "pheatmap"))
library(tidyverse)
library(pheatmap)

#Read the gene_counts.txt from featureCounts
# Read the table, skipping the comment lines that start with '#'
fc <- read.delim("counts/gene_counts.txt", comment.char = "#")

head(fc)[, 1:8]   # just take a peek at the first few columns




##make a clean count matrix:
# Keep gene IDs as row names
rownames(fc) <- fc$Geneid

# Keep only the count columns (samples)
countdata <- fc[, grep("SRR", colnames(fc))]

head(countdata)



#Rename the samples + create the sample info table
colnames(countdata)



#Rename:
colnames(countdata) <- c("CTRL1", "CTRL2", "siMYC1", "siMYC2")
colnames(countdata)




#create the metadata (colData) DESeq2 needs:
condition <- factor(c("CTRL", "CTRL", "siMYC", "siMYC"))

coldata <- data.frame(
  row.names = colnames(countdata),
  condition = condition
)

coldata


all(rownames(coldata) == colnames(countdata))


library(DESeq2)
BiocManager::install("DESeq2")
library(DESeq2)

DESeqDataSetFromMatrix


##Build the DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData   = coldata,
  design    = ~ condition
)
dds





#Filter out very low-count genes (optional but recommended)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
dds


#Run DESeq (normalisation + dispersion + testing)
dds <- DESeq(dds)



#Extract results: siMYC vs CTR
res <- results(dds, contrast = c("condition", "siMYC", "CTRL"))
summary(res)


#Sort by significance:
res_ordered <- res[order(res$padj), ]
head(res_ordered)




#Save the DE results to a file
#Create a folder
dir.create("DE_results", showWarnings = FALSE)





#Save
res_df <- as.data.frame(res_ordered)

write.csv(
  res_df,
  file = "DE_results/DESeq2_siMYC_vs_CTRL_results.csv"
)





#Quick QC plots 
#MA plot
plotMA(res, ylim = c(-5, 5))



#Variance stabilizing transform + PCA
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "condition")



#Heatmap of top 50 DE genes
##Install + load heatmap tools
install.packages("pheatmap")
install.packages("RColorBrewer")   # optional but nice

library(pheatmap)
library(RColorBrewer)


#Pick top 50 significant genes
# Remove genes with NA padj
res_clean <- res[!is.na(res$padj), ]

# Order by significance
res_ordered <- res_clean[order(res_clean$padj), ]

# Take top 50 DE genes
top50_genes <- rownames(res_ordered)[1:50]
top50_genes[1:5]  # just to check


library(DESeq2)

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# Remove genes with NA padj
res_clean <- res[!is.na(res$padj), ]

# Order by significance
res_ordered <- res_clean[order(res_clean$padj), ]

# Take top 50 DE genes
top50_genes <- rownames(res_ordered)[1:50]
top50_genes[1:5]  # quick check


vsd_mat <- assay(vsd)    # matrix of variance-stabilized counts



#Subset only top-50 genes
mat_top50 <- vsd_mat[top50_genes, ]


#Z-score rows for better heatmap visualization
mat_top50_z <- t(scale(t(mat_top50)))


#Make the heatmap
pheatmap(
  mat_top50_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 9, "RdBu")))(255),
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Top 50 Differentially Expressed Genes"
)


#Saving all plots 

#Create a folder for images

dir.create("plots", showWarnings = FALSE)


#Save MA plot
png("plots/MA_plot.png", width=1800, height=1600, res=200)
plotMA(res, ylim = c(-5, 5))
dev.off()

#Save PCA plot
png("plots/PCA_plot.png", width=1800, height=1600, res=200)
plotPCA(vsd, intgroup = "condition")
dev.off()


#Save dispersion plot
png("plots/Dispersion_plot.png", width=1800, height=1600, res=200)
plotDispEsts(dds)
dev.off()


#Save heatmap (top 50 DE genes)
png("plots/Heatmap_top50_DE_genes.png", width=2000, height=2400, res=200)

pheatmap(
  mat_top50_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 9, "RdBu")))(255),
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()


#Save Volcano Plot
png("plots/Volcano_plot.png", width=2000, height=1800, res=200)

with(res, plot(log2FoldChange, -log10(padj), pch=20,
               col=ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "red", "gray"),
               main="Volcano Plot"))

abline(v=c(-1,1), col="blue", lty=2)
abline(h=-log10(0.05), col="blue", lty=2)

dev.off()






