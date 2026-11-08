# SCENIC+ Workflow Overview

This folder contains the complete SCENIC+ pipeline for enhancer-driven gene regulatory network
(eGRN) inference and in silico perturbation analysis of gastric epithelial cell types.

The workflow is organized into the following scripts, to be run in order:

---

## Scripts

### `0_Subset_Epithelial_Clusters.ipynb`
Subsets epithelial cells from the integrated Seurat/Signac object prior to running SCENIC+.

- Clusters retained:
  - Parietal, Parietal_25, Foveolar, Dividing, Isthmus_5, Isthmus_12, Isthmus_20,
    Neck, Chief, Tuft, D_cells, Ghrl and Chga_Hdc
- Exports processed snRNA-seq object as `.h5ad` for use in SCENIC+

---

### `01_SCENICplus_setup_R.ipynb`
R-based setup script for preparing SCENIC+ inputs from the Seurat/Signac object.

- Exports the processed snRNA-seq Seurat object with cell type annotations to `.h5ad` format
- Prepares raw ATAC-seq fragment files for downstream pycisTopic preprocessing

---

### `02_SCENICplus_setup_Python.ipynb`
Python-based setup and eGRN inference script.

#### snATAC-seq preprocessing with pycisTopic
- Raw fragment files preprocessed using pycisTopic following author recommendations
- Topic modelling performed with MALLET, testing models with the following number of topics:
  `30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110`
- Final model selected based on model selection metrics: **80 topics**

#### Custom cistarget database
- Built using the `create_cisTarget_databases` Python package
- Inputs:
  - ATAC-seq consensus peaks calculated per cell type
  - Annotated motif collection provided by the SCENIC+ authors
  - Mouse TF list: `https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt`

#### eRegulon inference
- SCENIC+ run using default parameters
- A total of **468 eRegulons** were identified:
  - 170 direct
  - 298 extended

- eRegulon categories:
  - `+/+`: TF activates gene expression while opening chromatin
  - `-/+`: TF represses gene expression despite open chromatin at target genes

#### eRegulon filtering
- Filtering applied using `apply_std_filtering_to_eRegulons`:
  - Removal of negative region-to-gene (R2G) correlations (repressive regulatory relationships)
  - Extended eRegulons removed if a direct version exists (to prevent redundancy)
- After filtering: **160 eRegulons** retained (94 direct, 66 extended)
- Only positive (`+/+`) eRegulons used for downstream cell type-specific eGRN construction

#### eRegulon enrichment scoring
- AUC (Area Under the Curve) scores calculated per cell using the `AUCell` algorithm:
  - Gene-based AUC: enrichment of eRegulon target genes
  - Region-based AUC: enrichment of eRegulon target regions
- Regulon Specificity Score (RSS) calculated using Jensen-Shannon divergence to quantify
  lineage-specific activity of each eRegulon
- AUC and RSS values used jointly to link eRegulons to cell types

---

### `03_SCENICplus_Snakemake_config.yaml`
Snakemake configuration file for running the SCENIC+ pipeline.

- Defines parameters for the SCENIC+ snakemake pipeline run
- Run using default parameters

---

### `04_Perturbation-Model.ipynb`
In silico TF perturbation analysis using the SCENIC+ predictive modeling framework.

#### Perturbation approach
- For each gene in the eGRN, a regression model was trained to predict expression
  from upstream TFs with enhancer-linked binding sites
- TF knockdowns simulated by:
  1. Setting the expression of targeted TFs to zero
  2. Iteratively propagating predicted expression values through the network for **5 cycles**
     to capture indirect regulatory effects
- Simulated expression matrix co-embedded with the original SCENIC+ UMAP to visualize
  predicted cell-state shifts
- Delta-embedding used to quantify direction and magnitude of perturbation-induced transitions


## Key outputs

| Output | Description |
|--------|-------------|
| eRegulons | TF-to-enhancer-to-gene regulatory relationships |
| RSS heatmap | Cell type specificity of each eRegulon |
| AUC heatmap | Enrichment of eRegulon target genes and regions per cell |
| eGRN network plots | Top 10 edges per TF based on triplet ranking values |
| Perturbation UMAPs | Predicted cell-state shifts upon TF knockdown |
| Perturbation heatmaps | log2FC between original and perturbed gene expression |

---

## Key packages and resources

| Item | Details |
|------|---------|
| SCENIC+ | Bravo González-Blas et al., Nat Methods 2023 |
| pycisTopic | snATAC-seq preprocessing and topic modelling |
| MALLET | Topic modelling backend for pycisTopic |
| AUCell | eRegulon enrichment scoring |
| Mouse TF list | `https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt` |
| SCENIC+ snakemake pipeline | Default parameters |
