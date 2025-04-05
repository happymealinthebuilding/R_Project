# Results and Explanations

## 1. Expression Distribution Graphic

```
# the code for the graphic:
melted_df <- melt(airway_df[1:100, ])  # First 100 genes for simplicity

ggplot(melted_df, aes(x = variable, y = value)) +
  geom_boxplot(fill = "pink") +
  theme_minimal() +
  labs(title = "Expression Distribution", x = "Samples", y = "Expression")
```
![Expression Distribution](/Users/azratuncay/Desktop/R Projects/Expression Distribution.png)

### Explanation 
This block creates a boxplot showing the distribution of expression values for the first 100 genes across all samples.

## 2. Heat Map

```
# the code for the heat map:
var_genes <- apply(airway_df[, -1], 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:50]

heat_data <- as.matrix(airway_df[airway_df$gene %in% top_genes, -1])
rownames(heat_data) <- airway_df[airway_df$gene %in% top_genes, "gene"]

pheatmap(heat_data, scale = "row", clustering_distance_rows = "euclidean",
         clustering_method = "complete", show_rownames = FALSE)
```

![Heat Map](/Users/azratuncay/Desktop/R Projects/Heat Map.png)

### Explanation 
This section generates a heatmap of the top 50 most variable genes. The expression values are scaled by row and clustered using hierarchical clustering.

## 3. Volcano Plot: Treated vs Control

```
# the code for the volcano plot:
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = 'Volcano Plot: Treated vs Control')

```

![Volcano Plot: Treated vs Control](/Users/azratuncay/Desktop/R Projects/Volcano Plot: Treated vs Control.png)

### Explanation 
This block uses the EnhancedVolcano package to create a volcano plot. It visualizes the differentially expressed genes based on log2 fold change and p-value between the treated and control groups.

## 4.GO Enrichment - Biological Processes

```
# the code for the GO Enrichment Analysis:

ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP", # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

dotplot(ego, showCategory = 10)

```

![GO Enrichment- Biological Processes](/Users/azratuncay/Desktop/R Projects/GO Enrichment- Biological Processes.png)

### Explanation 
This section performs Gene Ontology (GO) enrichment analysis on the differentially expressed genes (DEGs), focusing on Biological Processes (BP). The dotplot function visualizes the top enriched biological processes.