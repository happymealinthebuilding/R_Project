# Basic Bioinformatics Analysis in R

This repository showcases a basic-level bioinformatics project focused on gene expression analysis using R. The analysis uses the "airway" dataset to explore differentially expressed genes, generate visualizations, and perform GO enrichment analysis. The project demonstrates various R techniques and packages, including DESeq2, ggplot2, pheatmap, and clusterProfiler.

## Installation

To run this project, you'll need to install the necessary R packages. The following script installs the required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("airway", "DESeq2", "ggplot2", "pheatmap", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"))
```

## Overview of the Project

### 1. Data Loading and Preprocessing
The airway dataset is loaded and transformed into a dataframe.
Gene names are added as a new column for better clarity.

### 2. Data Visualization
A boxplot is generated to visualize the expression distribution across samples using ggplot2.
A heatmap is created using the top 50 most variable genes to explore gene expression patterns.

### 3. Differential Expression Analysis
DESeq2 is used to perform differential expression analysis between treated and control groups.

A volcano plot is created using the EnhancedVolcano package to visualize differentially expressed genes.

### 4. Gene Ontology (GO) Enrichment
Gene symbols are mapped from Ensembl IDs.
GO enrichment analysis is performed using the clusterProfiler package, specifically focusing on Biological Processes.

## Results

The volcano plot highlights significant differentially expressed genes.
A GO enrichment bar plot and dot plot are generated to explore the enriched biological processes for the differentially expressed genes.

## Future Improvements

The project can be expanded to include more complex datasets.
Additional statistical tests and analyses could be integrated for a deeper exploration of gene expression data.
