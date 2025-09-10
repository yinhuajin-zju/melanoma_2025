# =============================================================================
# Method 3: Harmony Integration
# Author: Hanzhang
# Description: Integration of four Method 3 samples (0045, 0051, 0060, 0067)
#              using Harmony batch correction
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(harmony)
  library(future)
  library(stringr)
})

# =============================================================================
# Environment Configuration
# =============================================================================

set.seed(123)
options(future.globals.maxSize = 35 * 1024^3)
plan(multicore, workers = 4)

# =============================================================================
# Data Loading
# =============================================================================

M3_0045 <- readRDS("/path/to/M3_0045_wnn.RDS")
M3_0051 <- readRDS("/path/to/M3_0051_wnn.RDS")
M3_0060 <- readRDS("/path/to/M3_0060_wnn.RDS")
M3_0067 <- readRDS("/path/to/M3_0067_wnn.RDS")

seurat_list <- list(
  M3_0051 = RenameCells(M3_0051, add.cell.id = "M3_0051"),
  M3_0060 = RenameCells(M3_0060, add.cell.id = "M3_0060"),
  M3_0067 = RenameCells(M3_0067, add.cell.id = "M3_0067"),
  M3_0045 = RenameCells(M3_0045, add.cell.id = "M3_0045")
)

# =============================================================================
# Individual Sample Preprocessing
# =============================================================================

obj.list <- list(M3_0045, M3_0051, M3_0060, M3_0067)
names(obj.list) <- c("M3_0045", "M3_0051", "M3_0060", "M3_0067")

for (i in seq_along(obj.list)) {
  cat(paste("Processing", names(obj.list)[i], "...\n"))
  obj.list[[i]] <- SCTransform(
    obj.list[[i]],
    verbose = FALSE,
    new.assay.name = "SCT"
  )
}

# =============================================================================
# Sample Merging
# =============================================================================

merged <- merge(
  x = obj.list[[1]],
  y = obj.list[-1],
  add.cell.ids = names(obj.list),
  project = "M3_Harmony"
)

# =============================================================================
# Feature Selection and PCA
# =============================================================================

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
VariableFeatures(merged) <- features

merged <- RunPCA(
  merged,
  assay = "SCT",
  npcs = 50,
  verbose = FALSE
)

# =============================================================================
# Sample Identity Assignment
# =============================================================================

barcodes <- colnames(merged)
sample_label <- vapply(strsplit(barcodes, "_"),
                       function(x) paste(x[1:2], collapse = "_"),
                       FUN.VALUE = character(1))

merged$orig.ident <- sample_label

cat("Sample distribution:\n")
print(table(sample_label))

# =============================================================================
# Harmony Batch Correction
# =============================================================================

merged <- RunHarmony(
  object = merged,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:30,
  assay.use = "SCT",
  project.dim = FALSE,
  verbose = TRUE
)

# =============================================================================
# Post-Integration Analysis
# =============================================================================

merged <- RunUMAP(
  merged,
  reduction = "harmony",
  dims = 1:30
)

merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims = 1:30
)

merged <- FindClusters(
  merged,
  resolution = 0.5
)

# =============================================================================
# Visualization
# =============================================================================

cluster_plot <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
print(cluster_plot)

sample_plot <- DimPlot(merged, reduction = "umap", group.by = "orig.ident")
print(sample_plot)

# =============================================================================
# WNN Analysis and Save
# =============================================================================

seu <- RunWNNAnalysis(merged)
saveRDS(seu, file = "/path/to/outputs/M3_integrated.rds")


