# Pre-processing script to run on individual datasets

## Set up
params <- list(
  data = "/mnt/scratch/maithe/Gastric_Origins_Multiome/Project_GCG_15764_B01_SCR_Lane",
  sample_sequencing_id = "Pilot1_GCG_15614_B01_SCR_Lane",
  sample = "Pilot1",
  frags_path = "/mnt/scratch/maithe/Gastric_Origins_Multiome/Project_GCG_15614_B01_SCR_Lane/Sample_Ms-Stom-Pool-Gastric-Grant-Pilot-1-GEX/compbio/cellranger_override/outs/atac_fragments.tsv.gz",
  barcode_metrics_path = "/mnt/scratch/maithe/Gastric_Origins_Multiome/Project_GCG_15614_B01_SCR_Lane/Sample_Ms-Stom-Pool-Gastric-Grant-Pilot-1-GEX/compbio/cellranger_override/outs/per_barcode_metrics.csv",
  GEX_sample = "Sample_Ms-Stom-Pool-Gastric-Grant-Pilot-1-GEX",
  HTO_sample = "Sample_Ms-Stom-Pool-Gastric-Grant-Pilot-1-HTO",
  gex_dir = "/mnt/scratch/maithe/Gastric_Origins_Multiome/Project_GCG_15614_B01_SCR_Lane/Sample_Ms-Stom-Pool-Gastric-Grant-Pilot-1-GEX/compbio/cellranger_override/outs/filtered_feature_bc_matrix",
  hto_dir = "/mnt/scratch/maithe/Gastric_Origins_Multiome/Project_GCG_15614_B01_SCR_Lane/Sample_Ms-Stom-Pool-Gastric-Grant-Pilot-1-HTO/compbio/kallisto_bus_out/featurecounts",
  res_dir = "/mnt/data/maithe_2025/Multiome_Mouse_Healthy_2025/seurat/individual_preprocessing/results",
  table_dir = "tables",
  so_dir = "sobjs",
  fig_dir = "figures",
  metrics = "summary.csv",
  org = "Ms",
  macs3_path = "/home/maithe_2025/apps/miniconda3/envs/macs3/bin/macs3")
  
  
res_dir <- here(paste0(params$res_dir))
so_dir <- here(res_dir, params$so_dir)
table_dir <- here(res_dir, params$table_dir)
fig_dir <- here(res_dir, params$fig_dir)
hto_dir <- params$hto_dir # already determined in the params list
gex_dir <- params$gex_dir # already determined in the params list


packages <- c("Seurat","Signac", "here",
  "cowplot","tidyverse", "patchwork", "ggplot2", 
  "DoubletFinder",
  "clustree", "sceasy",
  "Matrix", "matrixStats",
  "Herper", "reticulate")

purrr::walk(packages, library, character.only = TRUE) 

## Genomic information
if (params$org == "Ms") { 
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
} else { # get human information 
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
}

seqlevelsStyle(annotations) <- "UCSC" # instead of "NCBI"

if (params$org == "Ms") { # for mouse 
  blacklists <- blacklist_mm10
  ribo_genes <- c("Gm42418", "AY036118")
  mt_str <- "^mt-"
  s.features <- tryCatch(convertHumanGeneList(Seurat::cc.genes.updated.2019$s.genes), error = function(e) {
    return(str_to_title(Seurat::cc.genes.updated.2019$s.genes))
  })
  g2m.features <- tryCatch(convertHumanGeneList(Seurat::cc.genes.updated.2019$g2m.genes), error = function(e) {
    return(str_to_title(Seurat::cc.genes.updated.2019$g2m.genes))
  })
} else { # if not working with mouse, use human info instead
  blacklists <- blacklist_hg38
  ribo_genes <- c("Gm42418", "AY036118")
  mt_str <- "^MT-"
  s.features <- Seurat::cc.genes.updated.2019$s.genes
  g2m.features <- Seurat::cc.genes.updated.2019$g2m.genes
}

## RNA
seurat_raw <- Read10X(data.dir = gex_dir)

seurat <- CreateSeuratObject(seurat_raw$`Gene Expression`, project = params$GEX_sample)
seurat$orig.ident <- params$GEX_sample
seurat <- RenameCells(seurat, new.names = str_remove(Cells(seurat), ".+_")) # renaming cells to keep just the barcodes 

### % MT genes 
DefaultAssay(seurat) <- "RNA"

seurat <- seurat %>% PercentageFeatureSet(
  pattern  = mt_str,
  col.name = "Percent_mito"
)

### % Ribo genes 
ribo_genes <- intersect(ribo_genes, rownames(seurat[["RNA"]]$counts))
if (length(ribo_genes) > 0) {
  seurat$Percent_ribo <- colSums(seurat[["RNA"]]$counts[ribo_genes, ]) / colSums(seurat[["RNA"]]$counts) * 100
}

### Cell cycle 
seurat <- seurat %>% CellCycleScoring(
  s.features = s.features,
  g2m.features = g2m.features,
  nbin = 15
)

## ATAC

frags <- CreateFragmentObject(
  path = params$frags_path,
  cells = Cells(seurat)
)

head(colnames(seurat_raw$Peaks))
colnames(seurat_raw$Peaks) <- str_remove(colnames(seurat_raw$Peaks), ".+_")

ATAC <- CreateChromatinAssay(
  counts = seurat_raw$Peaks,
  sep = c(":", "-"),
  fragments = frags,
  annotation = annotations
)

seurat[["ATAC"]] <- ATAC

### Read in cellranger metrics
metrics <- read_csv(params$barcode_metrics_path) %>%
  column_to_rownames("barcode") %>% # make barcode column as rownames 
  .[Cells(seurat), ] # subsets the metrics information based on cell barcodes 

max(metrics$is_cell) # check if at least one barcode was called a cell

### Find out the percentage of reads in peaks 
seurat$peak_region_fragments <- metrics$atac_peak_region_fragments
seurat$pct_reads_in_peaks <- metrics$atac_peak_region_fragments / metrics$atac_fragments * 100

### Blacklist ratio
seurat$blacklist_fraction <- FractionCountsInRegion(
  object = seurat,
  assay = "ATAC",
  regions = blacklists)
  
### Nucleosome signal and TSS enrichment 
DefaultAssay(seurat) <- "ATAC"
seurat <- NucleosomeSignal(seurat, assay = "ATAC") 
seurat <- TSSEnrichment(seurat, assay = "ATAC", fast = FALSE) 

### Recall peaks using MACS3
DefaultAssay(seurat) <- "ATAC"

peaks <- CallPeaks(seurat, 
                   macs2.path = params$macs3_path,
                   outdir = res_dir)

### Remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklists, invert = TRUE)

### Quantify counts in each peak
macs3_counts <- FeatureMatrix(
  fragments = frags,
  features = peaks,
  cells = colnames(seurat))
  
seurat[["peaks"]] <- CreateChromatinAssay(
  counts = macs3_counts,
  fragments = frags,
  annotation = annotations)
  
## QC

seurat@meta.data$orig.ident2 <- "Pilot1"

QC <- VlnPlot(
  object = seurat,
  features = c("nCount_RNA","nFeature_RNA", "Percent_mito", "Percent_ribo",
               "nCount_ATAC", "nFeature_ATAC", "nCount_peaks", "nFeature_peaks", 
               "peak_region_fragments", "TSS.enrichment", "nucleosome_signal", 
               "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "blacklist_fraction"), 
    group.by = "orig.ident2", pt.size = 0, ncol = 7) +  theme(axis.title.x = NULL)

ggsave(file.path(fig_dir, paste0(params$sample, "_QC_before_demux_and_filtering.png")), QC, width = 25, height = 15, units = "in")


qc_summary <- tibble(
  Metric = c(
    "nCount_RNA (min)", "nCount_RNA (max)",
    "nFeature_RNA (min)", "nFeature_RNA (max)",
    "Cells with <500 genes", "Cells with <100 genes", 
   "Cells with >8,000 genes", "Cells with >10,000 genes", "Cells with >15,000 genes",
    "Cells with >5% mito", "Cells with >10% mito", "Cells with >15% mito", "Cells with >20% mito",
      
    "nCount_ATAC (min) (cellranger)", 
    "nCount_peaks (min) (macs3)", 
      
    "Cells with <1000 peaks (cellranger)",
    "Cells with <1000 peaks (macs3)",

    "Cells with <2000 peaks (cellranger)", 
    "Cells with <2000 peaks (macs3)", 
      
    "nCount_ATAC (max) (cellranger)",
    "nCount_peaks (max) (macs3)",

    "Cells with >50,000 peaks (cellranger)", 
    "Cells with >50,000 peaks (macs3)",

    "Cells with >100,000 peaks (cellranger)",
    "Cells with >100,000 peaks (macs3)",
      
    "Cells with nucleosome_signal > 2",
    "Cells with TSS.enrichment < 1",
    "Cells with pct_reads_in_peaks < 15",
    "Cells with blacklist_fraction > 0.05"
  ),
  Value = c(
    min(seurat[["nCount_RNA"]]), max(seurat[["nCount_RNA"]]),
    min(seurat[["nFeature_RNA"]]), max(seurat[["nFeature_RNA"]]),
    sum(seurat[["nFeature_RNA"]] < 500), sum(seurat[["nFeature_RNA"]] < 100), 
    sum(seurat[["nFeature_RNA"]] > 8000), sum(seurat[["nFeature_RNA"]] > 10000), sum(seurat[["nFeature_RNA"]] > 15000),
    sum(seurat[["Percent_mito"]] > 5), sum(seurat[["Percent_mito"]] > 10), sum(seurat[["Percent_mito"]] > 15), sum(seurat[["Percent_mito"]] > 20),
      
    min(seurat[["nCount_ATAC"]]), 
    min(seurat[["nCount_peaks"]]), 

    sum(seurat[["nCount_ATAC"]] < 1000),
    sum(seurat[["nCount_peaks"]] < 1000),

    sum(seurat[["nCount_ATAC"]] < 2000), 
    sum(seurat[["nCount_peaks"]] < 2000), 
      
    max(seurat[["nCount_ATAC"]]),
    max(seurat[["nCount_peaks"]]),

    sum(seurat[["nCount_ATAC"]] > 50000), 
    sum(seurat[["nCount_peaks"]] > 50000),
      
    sum(seurat[["nCount_ATAC"]] > 100000),
    sum(seurat[["nCount_peaks"]] > 100000),
      
    sum(seurat[["nucleosome_signal"]] > 2),
    sum(seurat[["TSS.enrichment"]] < 1),
    sum(seurat[["pct_reads_in_peaks"]] < 15),
    sum(seurat[["blacklist_fraction"]] > 0.05)
  )
)

write.csv(qc_summary, row.names = F,
          here(table_dir, paste0(params$sample, "_QC_before_demux_and_filtering_summary.csv")))
          

## Filter

params_filtering <- list(
  gene_min = 100,  # min n of genes detected per cell                  
  gene_max = 20000, # max n of genes detected per cell                   
  mito_max = 10,  # max % of mito reads allowed per cell
  min_nCount_peaks = 1000, 
  max_nCount_peaks = 100000,
  nucleosome_signal = 2,
  tss_enrichment = 1,
  pct_reads_in_peaks = 15,
  blacklist_fraction = 0.05
)

seurat$rna_qc <- seurat@meta.data %>%
  mutate(rna_qc = case_when(
    Percent_mito > params_filtering$mito_max ~ "high_mito_reads",
    nFeature_RNA > params_filtering$gene_max ~ "high_gene_count",
    nFeature_RNA < params_filtering$gene_min ~ "low_gene_count",
    TRUE ~ "pass"
)) %>% pull(rna_qc)

seurat$atac_qc <- seurat@meta.data %>%
  mutate(atac_qc = case_when(
    nCount_peaks < params_filtering$min_nCount_peaks ~ "low_peak",
    nCount_peaks > params_filtering$max_nCount_peaks ~ "high_peak",
    pct_reads_in_peaks < params_filtering$pct_reads_in_peaks ~ "low_peak_pct",
    nucleosome_signal > params_filtering$nucleosome_signal ~ "high_nucleosome",
    TSS.enrichment < params_filtering$tss_enrichment ~ "low_TSS_enrichment",
    blacklist_fraction > blacklist_fraction ~ "high_blacklist",
    TRUE ~ "pass"
    )) %>% pull(atac_qc)
    
seurat_filt <- subset(seurat,
                      subset = rna_qc == "pass" & atac_qc == "pass")
                      
QC_filt <- VlnPlot(
  object = seurat_filt,
  features = c("nCount_RNA","nFeature_RNA", "Percent_mito", "Percent_ribo",  "S.Score", "G2M.Score", "nCount_ATAC", 
               "nFeature_ATAC", "nCount_peaks", "nFeature_peaks", "peak_region_fragments", "TSS.enrichment", 
               "nucleosome_signal", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "blacklist_fraction"),
  pt.size = 0, ncol = 7, group.by = "orig.ident2") + theme(axis.title.x = NULL)


ggsave(file.path(fig_dir, paste0(params$sample, "_QC_after_filtering_before_demux.png")), QC_filt, width = 25, height = 15, units = "in")

## Demultiplexing and Doublet removal

### Access hashing matrix 
hto_matrix <- readMM(here(hto_dir, "matrix.mtx"))
rownames(hto_matrix) <- readLines(here(hto_dir, "barcodes.tsv")) %>% str_c("-1")
colnames(hto_matrix) <- readLines(here(hto_dir, "features.tsv"))
colnames(hto_matrix) <- gsub("_", "-", colnames(hto_matrix))
hto_matrix <- as.matrix(hto_matrix)

hto_assay_filt <- CreateAssayObject(t(hto_matrix[Cells(seurat_filt),]), assay = "HTO_assay")
seurat_filt[["HTO"]] <- hto_assay_filt

### HTODemux LibNorm
hto_matrix_norm <- sweep(hto_matrix, 1, rowSums(hto_matrix), FUN = "/")
hto_matrix_norm <- round(hto_matrix_norm * 10000)
  
hto_assay_norm_filt <- CreateAssayObject(t(hto_matrix_norm[Cells(seurat_filt),]), assay = "HTO_LibNorm")
seurat_filt[["HTO_LibNorm"]] <- hto_assay_norm_filt

DefaultAssay(seurat_filt) <- "HTO_LibNorm"
seurat_filt <- NormalizeData(seurat_filt, assay = "HTO_LibNorm", normalization.method = "CLR")
seurat_filt <- HTODemux(seurat_filt, assay = "HTO_LibNorm", positive.quantile = 0.99, nsamples = 100)
knitr::kable(table(seurat_filt$HTO_LibNorm_classification.global))

Idents(seurat_filt) <- "HTO_LibNorm_classification.global" 

fl1 <- HTOHeatmap(seurat_filt, assay = "HTO_LibNorm", ncells = 5000) + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + labs(title = paste0(params$sample, "_LibNorm after filtering")) 
fl2 <- VlnPlot(seurat_filt, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

Idents(seurat_filt) <- "HTO_LibNorm_classification" 

fl3 <- VlnPlot(seurat_filt, 
        idents = c("tag-3", "tag-5", "tag-6", "tag-7", "tag-8", "tag9", "Negative"), 
        features = "nCount_RNA", pt.size = 0.1, log = TRUE)


ggsave(file.path(fig_dir, paste0(params$sample,"_LibNorm_after_filtering.png")), fl1+fl2+fl3, width = 30, height = 10, units = "in")
fl1+fl2+fl3 


### Inspecting CellCycle
if (params$org == "Ms") { # for mouse 
  s.features <- tryCatch(convertHumanGeneList(Seurat::cc.genes.updated.2019$s.genes), error = function(e) {
    return(str_to_title(Seurat::cc.genes.updated.2019$s.genes))
  })
  g2m.features <- tryCatch(convertHumanGeneList(Seurat::cc.genes.updated.2019$g2m.genes), error = function(e) {
    return(str_to_title(Seurat::cc.genes.updated.2019$g2m.genes))
  })
} else { # if not working with mouse, use human info instead
  s.features <- Seurat::cc.genes.updated.2019$s.genes
  g2m.features <- Seurat::cc.genes.updated.2019$g2m.genes
}

DefaultAssay(seurat_filt) <- "SCT"
seurat_filt <- seurat_filt %>% CellCycleScoring(
  s.features = s.features,
  g2m.features = g2m.features,
  nbin = 15
)

cc<- DimPlot(seurat_filt, group.by = "Phase", reduction = "umap") + ggtitle("Cell Cycle Phase - Pilot 1 \n after quick preprocessing for exploration")
cc

### Keep Singlets
Idents(seurat_filt) <- "HTO_LibNorm_classification.global" 
seurat_singlet <- subset(seurat_filt, idents = "Singlet")


## Processing
params2 <- list(
    n_dims_rna = 1:50,
    n_dims_atac = 2:50
)

## RNA Processing

DefaultAssay(seurat_singlet) <- "RNA"
seurat_singlet <- SCTransform(seurat_singlet, new.assay.name = "SCT")
seurat_singlet <- RunPCA(seurat_singlet, verbose = FALSE, assay = "SCT", npcs = 50)

e1 <- ElbowPlot(seurat_singlet, ndims = 50, reduction = "pca")
e1
ggsave(file.path(fig_dir, paste0(params$sample, "_ElbowPlotRNA_singlets.png")), e1, width = 30, height = 15, units = "in")

seurat_singlet <- FindNeighbors(seurat_singlet, reduction = "pca", k.param = 30, dims = params2$n_dims_rna)
seurat_singlet <- FindClusters(seurat_singlet, resolution = c(0.2, 0.5, 0.8, 1), algorithm = 1) # Louvain
seurat_singlet <- RunUMAP(seurat_singlet, dims = params2$n_dims_rna, slot = "data", reduction.key = "UMAP_", assay = "SCT", n.neighbors = 30L)  

a1 <- DimPlot(seurat_singlet, group.by = "SCT_snn_res.0.2", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets SCT - res.0.2", paste0(params$n_dims_rna, "PCs"))
a2 <- DimPlot(seurat_singlet, group.by = "SCT_snn_res.0.5", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets SCT - res.0.5", paste0(params$n_dims_rna, "PCs"))
a3 <- DimPlot(seurat_singlet, group.by = "SCT_snn_res.0.8", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets SCT - res.0.8", paste0(params$n_dims_rna, "PCs"))
a4 <- DimPlot(seurat_singlet, group.by = "SCT_snn_res.1", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets SCT - res.1", paste0(params$n_dims_rna, "PCs"))
ggsave(file.path(fig_dir, paste0(params$sample, "_SCT_ClusterResolutions.png")), a1+a2+a3+a4, width = 20, height = 15, units = "in")

## ATAC Processing

DefaultAssay(seurat_singlet) <- "peaks" # select assay where peaks were called with MACS3

seurat_singlet <- RunTFIDF(seurat_singlet, assay = "peaks", method = 1) 

seurat_singlet <-  FindTopFeatures(seurat_singlet, min.cutoff = "q0", assay = "peaks")
seurat_singlet <- RunSVD(seurat_singlet, assay = "peaks", n = 50, reduction.key = "LSI_", reduction.name = "lsi") 
DepthCor(seurat_singlet)

e2 <- ElbowPlot(seurat_singlet, ndims = 50, reduction = "lsi") + labs(title = "Elbow Plot - Peaks Assay \nfollowing FindTopFeatures + TFIDF + SVD")
ggsave(file.path(fig_dir, paste0(params$sample, "_ElbowPlotPeaks_singlets.png")), e2 + dp, width = 15, height = 7, units = "in")

seurat_singlet <- FindNeighbors(seurat_singlet, reduction = "lsi", k.param = 30, dims = params2$n_dims_atac)
seurat_singlet <- FindClusters(seurat_singlet, resolution = c(0.2, 0.5, 0.8, 1), algorithm = 1) # Louvain            

seurat_singlet <- RunUMAP(seurat_singlet, reduction = 'lsi', dims = params2$n_dims_atac, reduction.name = "umap.peaks", reduction.key = "peaksUMAP")

b1 <- DimPlot(seurat_singlet, group.by = "peaks_snn_res.0.2", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets Peaks - res.0.2")
b2 <- DimPlot(seurat_singlet, group.by = "peaks_snn_res.0.5", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets Peaks - res.0.5")
b3 <- DimPlot(seurat_singlet, group.by = "peaks_snn_res.0.8", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets Peaks - res.0.8")
b4 <- DimPlot(seurat_singlet, group.by = "peaks_snn_res.1", label = TRUE, pt.size = 0.5) + labs(title = "P1 Singlets Peaks - res.1")

ggsave(file.path(fig_dir, paste0(params$sample, "_PeaksAssay_ClusterResolutions.png")), b1+b2+b3+b4, width = 20, height = 15, units = "in")

### Score Gene Activity
gene.activities.peaks <- GeneActivity(seurat_singlet, assay = "peaks", extend.upstream = 2000, extend.downstream = 0, biotypes = "protein_coding")

seurat_singlet[["GApeaks"]] <- CreateAssayObject(counts = gene.activities.peaks)

seurat_singlet <- NormalizeData(
  object = seurat_singlet,
  assay = "GApeaks",
  scale.factor = median(seurat_singlet$nCount_GApeaks))
  

## WNN Integration of RNA and ATAC 

DefaultAssay(seurat_singlet) <- "SCT"

multiome_singlet <- FindMultiModalNeighbors(
  object = seurat_singlet,
  reduction.list = list("pca", "lsi"),
  dims.list = list(params2$n_dims_rna, params2$n_dims_atac),
  verbose = TRUE, k.nn = 30)
  
multiome_singlet <- RunUMAP(
  object = multiome_singlet,
  nn.name = "weighted.nn",
  verbose = TRUE, 
  reduction.name = "jointumap",
  reduction.key = "jointUMAP_",
  n.neighbors = 30L
)

multiome_singlet <- FindClusters(
  multiome_singlet, 
  graph.name = "wsnn", 
  algorithm = 3, verbose = FALSE, resolution = c(0.2, 0.5, 0.8, 1))
  
Reductions(multiome_singlet)

p1 <- DimPlot(multiome_singlet, reduction = 'jointumap', group.by = "wsnn_res.0.2", label = TRUE, pt.size = 0.5, shuffle = T) + labs(title = "P1 WNN knn 30 - res.0.2")
p2 <- DimPlot(multiome_singlet, reduction = 'jointumap', group.by = "wsnn_res.0.5", label = TRUE, pt.size = 0.5, shuffle = T) + labs(title = "P1 WNN knn 30 - res.0.5")
p3 <- DimPlot(multiome_singlet, reduction = 'jointumap', group.by = "wsnn_res.0.8", label = TRUE, pt.size = 0.5, shuffle = T) + labs(title = "P1 WNN knn 30- res.0.8")
p4 <- DimPlot(multiome_singlet, reduction = 'jointumap', group.by = "wsnn_res.1", label = TRUE, pt.size = 0.5, shuffle = T) + labs(title = "P1 WNN knn 30- res.1")
p1 + p2 + p3 + p4

ggsave(file.path(fig_dir, paste0(params$sample, "_WNN_knn30_ClusterResolution.png")), p1+p2+p3+p4, width = 20, height = 15, units = "in")

### Clustree to inspect clustering resolution 
c1 <- clustree(multiome_singlet, prefix = paste0("SCT_snn_res."), prop_filter = 0.3)
c2 <- clustree(multiome_singlet, prefix = paste0("peaks_snn_res."), prop_filter = 0.3)
c3 <- clustree(multiome_singlet, prefix = paste0("wsnn_res."), prop_filter = 0.3)

ggsave(file.path(fig_dir, paste0(params$sample, "_Clustree_All_Assays.png")), c1 +c2 +c3, width = 30, height = 10, units = "in")

saveRDS(
  multiome_singlet,
  here(so_dir, paste0(params$sample, "_seurat_singlets_WNNmultiome.rds")))





