if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("airway", "DESeq2", "ggplot2", "pheatmap", "EnhancedVolcano"))
library(airway)
data("airway")
airway_df <- as.data.frame(assay(airway))
head(airway_df)
# Convert rownames to a column
airway_df$gene <- rownames(airway_df)
airway_df <- airway_df[, c(ncol(airway_df), 1:(ncol(airway_df)-1))]
summary(airway_df[, -1])
install.packages("reshape2")

library(ggplot2)

# Melt the data for ggplot
library(reshape2)
melted_df <- melt(airway_df[1:100, ])  # First 100 genes for simplicity

ggplot(melted_df, aes(x = variable, y = value)) +
  geom_boxplot(fill = "pink") +
  theme_minimal() +
  labs(title = "Expression Distribution", x = "Samples", y = "Expression")
library(pheatmap)

# Use top 50 most variable genes
var_genes <- apply(airway_df[, -1], 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:50]

heat_data <- as.matrix(airway_df[airway_df$gene %in% top_genes, -1])
rownames(heat_data) <- airway_df[airway_df$gene %in% top_genes, "gene"]

pheatmap(heat_data, scale = "row", clustering_distance_rows = "euclidean",
         clustering_method = "complete", show_rownames = FALSE)
library(DESeq2)

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds)

# View top differentially expressed genes
head(res[order(res$pvalue), ])
library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = 'Volcano Plot: Treated vs Control')
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"))
head(deg_genes)
keys(org.Hs.eg.db, keytype = "SYMBOL")
entrez_ids <- bitr(deg_genes,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
library(org.Hs.eg.db)

# Map Ensembl IDs to gene symbols
deg_genes_symbols <- mapIds(org.Hs.eg.db, keys = deg_genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# View the mapped gene symbols
head(deg_genes_symbols)

ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP", # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory = 10, title = "GO Enrichment - Biological Processes")

# OR try:
dotplot(ego, showCategory = 10)
install.packages("clusterProfiler")




