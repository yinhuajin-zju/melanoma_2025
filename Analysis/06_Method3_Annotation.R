# =============================================================================
# Method 3 integrated data annotation by using Celltypist
# Author: Hanzhang
# Description: Use Celltypist to do the automatic annotation and manually adjust
#              by using CellScoreModule with markers
# =============================================================================
library(Seurat)
library(stringr)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)

# =============================================================================
# Environment Setup
# =============================================================================

required_pkg <- c("Seurat")

for (pkg in required_pkg) {
  library(pkg, character.only = TRUE)
}

rm(pkg)

# =============================================================================
# Data Loading and Initial Exploration
# =============================================================================

seu <- readRDS("./Data/M3_integrated.rds")

head(seu@meta.data)
DimPlot(seu, group.by = "seurat_clusters", reduction = "umap")
DimPlot(seu, group.by = "orig.ident", reduction = "umap")

DefaultAssay(seu) <- "RNA"
table(rownames(seu))

barcodes <- colnames(seu)
sample_label <- str_extract(barcodes, "^[^_]+_[^_]+")
seu$orig.ident <- sample_label

# =============================================================================
# Expression Matrix Export
# =============================================================================

seu <- JoinLayers(seu, assay = "RNA")
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")

expr_mat <- t(counts)
expr_df <- as.data.frame(as.matrix(expr_mat))
expr_df <- cbind(cell = rownames(expr_df), expr_df)

fwrite(expr_df,
       file = "seurat_raw_counts_cells_rows.csv",
       row.names = FALSE,
       quote = FALSE,
       sep = ",")

expr_mat2 <- counts
expr_df2 <- as.data.frame(as.matrix(expr_mat2))
expr_df2 <- cbind(gene = rownames(expr_df2), expr_df2)

fwrite(expr_df2,
       file = "seurat_raw_counts_genes_rows.csv.gz",
       row.names = FALSE,
       quote = FALSE,
       sep = ",")

# =============================================================================
# Cell Type Annotation Integration
# =============================================================================

annotations <- read.csv("predicted_labels.csv")
head(annotations)

colnames(annotations)[1] <- "cell_barcode"
rownames(annotations) <- annotations$cell_barcode

annotations_to_add <- annotations %>%
  select(predicted_labels, over_clustering, majority_voting)

head(annotations_to_add)

seu <- AddMetaData(
  object = seu,
  metadata = annotations_to_add
)

seu$majority_voting[seu$majority_voting == "Melanocyte"] <- "Melanoma cell"

head(seu@meta.data)

# =============================================================================
# Module Score-Based Cell Type Refinement
# =============================================================================

marker_list <- list(
  NK = c("KLRD1", "NCAM1", "KLRK1", "FGFBP2", "NCR1", "PRF1", "GZMB", "GNLY", "KLRB1", "NKG7"),
  Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD20")
)

seu$cell_anno_original <- seu$cell_anno

seu <- AddModuleScore(
  object = seu,
  features = marker_list,
  name = names(marker_list)
)

score_cols <- paste0(names(marker_list), 1:length(marker_list))

scores_df <- seu@meta.data[, score_cols]

max_score_cell_type <- apply(scores_df, 1, function(row) {
  if (max(row) > 0) {
    return(names(which.max(row)))
  } else {
    return(NA)
  }
})

max_score_cell_type <- gsub("[0-9]+$", "", max_score_cell_type)
cells_to_update <- !is.na(max_score_cell_type)
seu$cell_anno[cells_to_update] <- max_score_cell_type[cells_to_update]

# =============================================================================
# Cell Type Identity Setting and Visualization
# =============================================================================

Idents(seu) <- "majority_voting"
DimPlot(seu, reduction = "umap")
table(Idents(seu))

DimPlot(seu, reduction = "umap", label = TRUE)
table(Idents(seu))

FeaturePlot(seu, features = score_cols, label = TRUE)

# =============================================================================
# Manual Cell Type Corrections
# =============================================================================

seu$cell_anno[seu$cell_anno == "DC1"] <- "Melanoma cell"
seu$cell_anno[seu$cell_anno == "Macro_1"] <- "Macrophage"
seu$cell_anno[seu$cell_anno == "Bcell"] <- "B cell"
seu$cell_anno[seu$cell_anno == "Pericyte_1"] <- "Pericyte"
seu$cell_anno[seu$cell_anno == "VE1"] <- "Endothelial cell"
seu$cell_anno[seu$cell_anno == "Inf_mac"] <- "Inflammatory Macrophage"

DimPlot(seu, group.by = "cell_anno", label = TRUE)

marker_list <- list(
  NK = c("KLRD1", "NCAM1", "KLRK1", "FGFBP2", "NCR1", "PRF1", "GZMB", "GNLY", "KLRB1", "NKG7"),
  Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD20")
)

seu$cell_anno_original <- seu$cell_anno

seu <- AddModuleScore(
  object = seu,
  features = marker_list,
  name = names(marker_list)
)

score_cols <- paste0(names(marker_list), 1:length(marker_list))

print(score_cols)

if (!all(score_cols %in% colnames(seu@meta.data))) {
  stop("Generated score_cols do not match column names in meta.data!")
}

scores_df <- seu@meta.data[, score_cols]

max_score_cell_type <- apply(scores_df, 1, function(row) {
  if (max(row) > 0) {
    return(names(which.max(row)))
  } else {
    return(NA)
  }
})

max_score_cell_type <- gsub("[0-9]+$", "", max_score_cell_type)

head(max_score_cell_type)
table(max_score_cell_type, useNA = "ifany")

cells_to_update <- !is.na(max_score_cell_type)
seu$cell_anno[cells_to_update] <- max_score_cell_type[cells_to_update]

cat("\nOriginal annotation statistics:\n")
print(table(seu$cell_anno_original))

cat("\nUpdated annotation statistics:\n")
print(table(seu$cell_anno))

FeaturePlot(seu, features = score_cols, label = TRUE)
DimPlot(seu, group.by = "cell_anno", label = TRUE)

# =============================================================================
# Save Results
# =============================================================================

saveRDS(seu, "./Data/M3_integrated_v2_with_anno.rds")
