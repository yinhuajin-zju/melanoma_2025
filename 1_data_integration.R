library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(future)
library(stringr)

# =============================================================================
# Configuration & Setup
# =============================================================================

# Set reproducible seed and library path
set.seed(123) 
.libPaths("/home/pj/R/x86_64-pc-linux-gnu-library/4.3")

# Configure parallel processing
options(future.globals.maxSize = 35 * 1024^3)
plan(multicore, workers = 4)

# =============================================================================
# Data Loading & Preparation
# =============================================================================

# Define sample metadata - single source of truth
SAMPLES <- c("M3_0045", "M3_0051", "M3_0060", "M3_0067")
SAMPLE_FILES <- paste0("updated_", gsub("M3_", "", SAMPLES), "_wnn.RDS")

# Load all samples with consistent naming
seurat_objects <- setNames(
  lapply(SAMPLE_FILES, readRDS),
  SAMPLES
)

# Add cell identifiers - prevents cell name collisions
seurat_objects <- Map(
  function(obj, sample_id) RenameCells(obj, add.cell.id = sample_id),
  seurat_objects,
  SAMPLES
)

# =============================================================================
# Preprocessing Pipeline
# =============================================================================

# Apply SCTransform normalization to all samples
# Note: Skip if objects are already SCTransformed
seurat_objects <- lapply(seurat_objects, function(obj) {
  SCTransform(obj, verbose = FALSE, new.assay.name = "SCT")
})

# Merge all objects into single dataset
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_objects)

# Extract sample identity from cell barcodes
# Data flow: barcode -> sample_id -> metadata
cell_barcodes <- colnames(merged_seurat)
sample_labels <- vapply(
  strsplit(cell_barcodes, "_"),
  function(parts) paste(parts[1:2], collapse = "_"),
  FUN.VALUE = character(1)
)
merged_seurat$orig.ident <- sample_labels

# =============================================================================
# Integration Workflow
# =============================================================================

# Select integration features across all samples
integration_features <- SelectIntegrationFeatures(
  object.list = seurat_objects, 
  nfeatures = 3000
)
VariableFeatures(merged_seurat) <- integration_features

# Dimensionality reduction pipeline
merged_seurat <- merged_seurat %>%
  RunPCA(assay = "SCT", npcs = 50, verbose = FALSE) %>%
  RunHarmony(
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    dims.use = 1:30,
    assay.use = "SCT",
    verbose = TRUE
  ) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

# =============================================================================
# Quality Control Visualization
# =============================================================================

# Create diagnostic plots
cluster_plot <- DimPlot(
  merged_seurat, 
  reduction = "umap", 
  group.by = "seurat_clusters", 
  label = TRUE,
  pt.size = 0.5
) + labs(title = "Clusters")

batch_plot <- DimPlot(
  merged_seurat, 
  reduction = "umap", 
  group.by = "orig.ident",
  pt.size = 0.5
) + labs(title = "Batch Integration")

# Display integration quality
print(cluster_plot | batch_plot)

# =============================================================================
# Final Processing & Export
# =============================================================================

# Add WNN analysis for multimodal integration
integrated_seurat <- RunWNNAnalysis(merged_seurat)

# Save final integrated object
saveRDS(integrated_seurat, file = "./Data/M3_integrated_v2.rds")
