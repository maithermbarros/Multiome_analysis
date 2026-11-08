# Seurat/Signac Workflow Overview

This folder contains the complete Seurat/Signac-based preprocessing and integration pipeline for single-nucleus multiome (snRNA-seq + snATAC-seq) datasets.

The workflow is organized into two key stages, each with its own scripts:

---

## 1 - `seurat_preprocessing/`

Each individual dataset (Batch 1: 6 samples; Batch 2: 4 samples) is processed independently for both RNA and ATAC modalities.

### snRNA-seq preprocessing
- Quality control:
  - Removal of cells expressing fewer than 100 genes
  - Removal of cells with mitochondrial gene content exceeding 10%
- Normalization using `SCTransform`
- Dimensionality reduction via PCA (`RunPCA`, first 50 PCs)
- Graph-based clustering using the Louvain algorithm (`FindNeighbors`, `FindClusters`)
- UMAP embedding (`RunUMAP`)

### snATAC-seq preprocessing
- Quality control:
  - Removal of cells with fewer than 1,000 or more than 100,000 fragments
  - Removal of cells with transcription start site enrichment < 1
  - Removal of cells with percentage of reads in peaks < 15%
- Peak calling with MACS3 (via Signac), excluding genomic blacklist regions
- Normalization using Term Frequency–Inverse Document Frequency (`RunTFIDF`)
- Feature selection with `FindTopFeatures`
- Dimensionality reduction using Singular Value Decomposition (`RunSVD`, 50 LSI components)
- Graph-based clustering (`FindNeighbors`, `FindClusters`)
- UMAP embedding (`RunUMAP`)

### Sample demultiplexing
- HTO libraries processed using the Kallisto/Kite workflow
- Singlet identification using a modified `HTODemux` with hashtag library-size normalization

### Per-batch joint embedding
- Weighted Nearest Neighbors (`FindMultiModalNeighbors`) applied to jointly cluster RNA and ATAC modalities within each batch

---

## 2 - `seurat_integration/`

Once all datasets are preprocessed, they are integrated into a unified multimodal object.

### RNA integration
- RNA assays from both batches merged with the `merge` function
- Re-processed with `SCTransform`, `RunPCA`, `FindNeighbors`, `FindClusters`, `RunUMAP`
- Batch correction using Harmony (`IntegrateLayers`, method = `HarmonyIntegration`)

### ATAC integration
- Unified peak set created across both ATAC assays using `reduce` (GenomicRanges) and quantified with `FeatureMatrix`
- ATAC assays merged and re-processed (`FindTopFeatures`, `RunTFIDF`, `RunSVD`, `FindNeighbors`, `FindClusters`, `RunUMAP`)
- Integration anchors identified with `FindIntegrationAnchors` and embeddings integrated with `IntegrateEmbeddings`

### Joint multimodal clustering
- WNN applied to integrated RNA and ATAC assays (`FindMultiModalNeighbors`)
- UMAP embedding and clustering at four resolutions: `0.2`, `0.5`, `0.8`, and `1.0`

### Differential expression and cell annotation
- Pseudo-bulk aggregate expression computed per animal per cluster (`AggregateExpression`, assay = `"SCT"`)
- Differentially expressed genes identified using DESeq2 (`FindAllMarkers`, `test.use = "DESeq2"`)
- Manual cell annotation based on known marker genes and top DEGs per cluster

---

## Reference genome and key packages

| Item | Version/Details |
|------|----------------|
| Reference genome | GRCm38 (mm10) — cellranger-arc-mm10-2020-A |
| Seurat | v5.3.0 |
| Signac | v1.14.0 |
