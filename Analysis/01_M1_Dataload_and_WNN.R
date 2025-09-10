# =============================================================================
# Method 1 Analysis: Comparison of Fixed vs Non-fixed Samples
# Author: Hanzhang
# Description: Data load and analysis of scRNA-seq and scATAC-seq data from 
#              fixed and non-fixed Method 1 samples using Seurat and Signac
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(dittoSeq)
  library(patchwork)
  library(plot1cell)
  library(RColorBrewer)
})

# =============================================================================
# Configuration and File Paths
# =============================================================================

setwd("/path/to/your/project/directory")

data_paths <- list(
  nofix = "/path/to/data/S-MEL_GE_nofix/outs/filtered_feature_bc_matrix.h5",
  fix = "/path/to/data/S-MEL_GE_fix/outs/filtered_feature_bc_matrix.h5"
)

fragment_paths <- list(
  nofix = "/path/to/data/S-MEL_GE_nofix/outs/atac_fragments.tsv.gz",
  fix = "/path/to/data/S-MEL_GE_fix/outs/atac_fragments.tsv.gz"
)

# =============================================================================
# Data Import and Initial Processing
# =============================================================================

inputdata.nofix <- Read10X_h5(data_paths$nofix)
inputdata.fix <- Read10X_h5(data_paths$fix)

rna_counts_nofix <- inputdata.nofix$`Gene Expression`
atac_counts_nofix <- inputdata.nofix$Peaks

rna_counts_fix <- inputdata.fix$`Gene Expression`
atac_counts_fix <- inputdata.fix$Peaks

# =============================================================================
# Create Seurat Objects and Add Metadata
# =============================================================================

M1_fix <- CreateSeuratObject(counts = rna_counts_fix, project = "M1_fix")
M1_nofix <- CreateSeuratObject(counts = rna_counts_nofix, project = "M1_nofix")

M1_fix[["percent.mt"]] <- PercentageFeatureSet(M1_fix, pattern = "^MT-")
M1_nofix[["percent.mt"]] <- PercentageFeatureSet(M1_nofix, pattern = "^MT-")

M1_fix$orig.id <- "M1_fix"
M1_nofix$orig.id <- "M1_nofix"

cat("M1_fix metadata:\n")
head(M1_fix@meta.data)

cat("M1_nofix metadata:\n")
head(M1_nofix@meta.data)

# =============================================================================
# Quality Control Visualization (RNA only)
# =============================================================================

qc_features <- c("percent.mt", "nCount_RNA", "nFeature_RNA")

p_fix_rna <- VlnPlot(M1_fix, features = qc_features, pt.size = 0, cols = "red") +
  ggtitle("M1_fix - RNA QC Metrics")

p_nofix_rna <- VlnPlot(M1_nofix, features = qc_features, pt.size = 0, cols = "red") +
  ggtitle("M1_nofix - RNA QC Metrics")

print(p_fix_rna)
print(p_nofix_rna)

# =============================================================================
# ATAC Data Processing
# =============================================================================

process_atac_counts <- function(atac_counts) {
  
  grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  atac_counts_filtered <- atac_counts[as.vector(grange_use), ]
  
  return(atac_counts_filtered)
}

atac_counts_nofix <- process_atac_counts(atac_counts_nofix)
atac_counts_fix <- process_atac_counts(atac_counts_fix)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# =============================================================================
# Create Chromatin Assays and Add to Seurat Objects
# =============================================================================

create_chromatin_assay <- function(counts, fragment_file, annotations) {
  CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = fragment_file,
    min.cells = 10,
    annotation = annotations
  )
}

M1_fix[["ATAC"]] <- create_chromatin_assay(atac_counts_fix, fragment_paths$fix, annotations)
M1_nofix[["ATAC"]] <- create_chromatin_assay(atac_counts_nofix, fragment_paths$nofix, annotations)

saveRDS(M1_fix, "outputs/M1_fix_raw.RDS")
saveRDS(M1_nofix, "outputs/M1_nofix_raw.RDS")

# =============================================================================
# Comprehensive Quality Control Visualization
# =============================================================================

multiome_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC")

plot_fix_multiome <- VlnPlot(
  M1_fix, 
  features = multiome_features, 
  ncol = 5,
  log = TRUE, 
  pt.size = 0
) + 
  NoLegend() + 
  ggtitle("M1_fix - Multiome QC Metrics")

plot_nofix_multiome <- VlnPlot(
  M1_nofix, 
  features = multiome_features, 
  ncol = 5,
  log = TRUE, 
  pt.size = 0
) + 
  NoLegend() + 
  ggtitle("M1_nofix - Multiome QC Metrics")

print(plot_fix_multiome)
print(plot_nofix_multiome)

# =============================================================================
# Quality Control Filtering
# =============================================================================

filtering_criteria <- list(
  nCount_ATAC = c(5e3, 7e4),
  nCount_RNA = c(1000, 22000),
  percent.mt = c(0, 20)
)

M1_fix <- subset(
  x = M1_fix,
  subset = nCount_ATAC > filtering_criteria$nCount_ATAC[1] &
    nCount_ATAC < filtering_criteria$nCount_ATAC[2] &
    nCount_RNA > filtering_criteria$nCount_RNA[1] &
    nCount_RNA < filtering_criteria$nCount_RNA[2] &
    percent.mt < filtering_criteria$percent.mt[2]
)

M1_nofix <- subset(
  x = M1_nofix,
  subset = nCount_ATAC > filtering_criteria$nCount_ATAC[1] &
    nCount_ATAC < filtering_criteria$nCount_ATAC[2] &
    nCount_RNA > filtering_criteria$nCount_RNA[1] &
    nCount_RNA < filtering_criteria$nCount_RNA[2] &
    percent.mt < filtering_criteria$percent.mt[2]
)

cat(paste("After filtering - M1_fix:", ncol(M1_fix), "cells\n"))
cat(paste("After filtering - M1_nofix:", ncol(M1_nofix), "cells\n"))

# =============================================================================
# RNA Data Processing and Dimension Reduction
# =============================================================================

process_rna_data <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
  return(seurat_obj)
}

cat("Processing RNA data...\n")
M1_fix <- process_rna_data(M1_fix)
M1_nofix <- process_rna_data(M1_nofix)

# =============================================================================
# ATAC Data Processing and Dimension Reduction
# =============================================================================

process_atac_data <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "ATAC"
  seurat_obj <- RunTFIDF(seurat_obj) %>%
    FindTopFeatures(min.cutoff = "q0") %>%
    RunSVD() %>%
    RunUMAP(reduction = "lsi", dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  return(seurat_obj)
}

M1_fix <- process_atac_data(M1_fix)
M1_nofix <- process_atac_data(M1_nofix)

# =============================================================================
# Weighted Nearest Neighbor (WNN) Analysis
# =============================================================================

perform_wnn_analysis <- function(seurat_obj, k_nn = 20) {
  seurat_obj <- FindMultiModalNeighbors(
    seurat_obj, 
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:20, 2:20),
    k.nn = ifelse(k_nn == 20, 20, 10),  # Special parameter for M1_fix
    knn.range = ifelse(k_nn == 20, 200, 100)
  ) %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
    FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  
  return(seurat_obj)
}

cat("Performing WNN analysis...\n")
M1_nofix <- perform_wnn_analysis(M1_nofix, k_nn = 20)
M1_fix <- perform_wnn_analysis(M1_fix, k_nn = 10)  # Different parameters for fix sample

# =============================================================================
# Visualization: Basic UMAP Plots
# =============================================================================

create_basic_umap_plots <- function(seurat_obj, title_prefix) {
  p1 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste(title_prefix, "- RNA"))
  
  p2 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste(title_prefix, "- ATAC"))
  
  p3 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste(title_prefix, "- WNN"))
  
  return(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
}

basic_plots_nofix <- create_basic_umap_plots(M1_nofix, "M1_nofix")
basic_plots_fix <- create_basic_umap_plots(M1_fix, "M1_fix")

print(basic_plots_nofix)
print(basic_plots_fix)

# =============================================================================
# Differential Expression Analysis
# =============================================================================

find_all_markers_safe <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj[["RNA"]]$data <- as(object = seurat_obj[["RNA"]]$counts, Class = "dgCMatrix")
  
  all_markers <- FindAllMarkers(object = seurat_obj, verbose = TRUE)
  
  return(all_markers)
}


all_markers_M1_fix <- find_all_markers_safe(M1_fix)
all_markers_M1_nofix <- find_all_markers_safe(M1_nofix)

cat("Top markers for M1_fix:\n")
head(all_markers_M1_fix)

cat("Top markers for M1_nofix:\n")
head(all_markers_M1_nofix)

# =============================================================================
# Save Final Results
# =============================================================================


saveRDS(M1_fix, "outputs/M1_fix_wnn.RDS")
saveRDS(M1_nofix, "outputs/M1_nofix_wnn.RDS")


