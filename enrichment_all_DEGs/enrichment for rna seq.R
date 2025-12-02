setwd("C:/Users/megha/projects/MYC_Project/LAYER2_RNAseq")
getwd() 


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", 
                       "org.Hs.eg.db",
                       "enrichplot",

                                              "ReactomePA"))
BiocManager::install("reactome.db")
BiocManager::install("ReactomePA")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)


#Reload your DESeq2 results
res_df <- read.csv("DE_results/DESeq2_siMYC_vs_CTRL_results.csv", row.names = 1)
head(res_df)
#Extract ALL significant DEGs
deg_all <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
nrow(deg_all)



#Convert ENSG → SYMBOL → ENTREZ
clean_ids_all <- gsub("\\..*", "", rownames(deg_all))

symbols_all <- mapIds(org.Hs.eg.db,
                      keys = clean_ids_all,
                      keytype = "ENSEMBL",
                      column = "SYMBOL",
                      multiVals = "first")

symbols_all <- symbols_all[!is.na(symbols_all)]

entrez_all <- mapIds(org.Hs.eg.db,
                     keys = symbols_all,
                     keytype = "SYMBOL",
                     column = "ENTREZID",
                     multiVals = "first")

entrez_all <- entrez_all[!is.na(entrez_all)]
length(entrez_all)




#GO BP enrichment
ego_all <- enrichGO(
  gene         = entrez_all,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)




#KEGG enrichment
ekegg_all <- enrichKEGG(
  gene = entrez_all,
  organism = "hsa",
  pvalueCutoff = 0.05
)




#Reactome enrichment
ereact_all <- enrichPathway(
  gene = entrez_all,
  pvalueCutoff = 0.05,
  readable = TRUE
)






head(ego_all)
head(ekegg_all)
head(ereact_all)




#Create an “enrichment_all_DEGs” folder
dir.create("enrichment_all_DEGs", showWarnings = FALSE)


#Save GO results to CSV
write.csv(as.data.frame(ego_all),
          "enrichment_all_DEGs/GO_BP_all_DEGs.csv")


#Save KEGG results
write.csv(as.data.frame(ekegg_all),
          "enrichment_all_DEGs/KEGG_all_DEGs.csv")



#Save Reactome results
write.csv(as.data.frame(ereact_all),
          "enrichment_all_DEGs/Reactome_all_DEGs.csv")



#Make folder for images:
dir.create("enrichment_all_DEGs/plots", showWarnings = FALSE)

#GO dotplot
png("enrichment_all_DEGs/plots/GO_BP_dotplot.png", width=2000, height=1800, res=200)
dotplot(ego_all, showCategory = 20)
dev.off()


#KEGG dotplot
png("enrichment_all_DEGs/plots/KEGG_dotplot.png", width=2000, height=1800, res=200)
dotplot(ekegg_all, showCategory = 20)
dev.off()


#Reactome dotplot
png("enrichment_all_DEGs/plots/Reactome_dotplot.png", width=2000, height=1800, res=200)
dotplot(ereact_all, showCategory = 20)
dev.off()


#Save R objects (recommended)
saveRDS(ego_all, "enrichment_all_DEGs/GO_BP_all_DEGs.rds")
saveRDS(ekegg_all, "enrichment_all_DEGs/KEGG_all_DEGs.rds")
saveRDS(ereact_all, "enrichment_all_DEGs/Reactome_all_DEGs.rds")


