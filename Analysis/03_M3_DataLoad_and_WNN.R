# =============================================================================
# Method 3 Analysis: M3 data load and preprocess
# Author: Hanzhang
# Date: [Current Date]
# Description: Comprehensive analysis of scRNA-seq and scATAC-seq data from 
#              four different samples (0045, 0051, 0060, 0067) using Seurat 
#              and Signac workflows
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(dittoSeq)
  library(ggrastr)
})

library(Nebulosa)

# =============================================================================
# Configuration and File Paths
# =============================================================================

setwd("/path/to/your/method3/directory")

sample_info <- list(
  samples = c("0045", "0051", "0060", "0067"),
  data_paths = list(
    "0045" = "/path/to/data/Method3_0045/outs/filtered_feature_bc_matrix.h5",
    "0051" = "/path/to/data/Method3_0051/outs/filtered_feature_bc_matrix.h5",
    "0060" = "/path/to/data/Method3_0060/outs/filtered_feature_bc_matrix.h5",
    "0067" = "/path/to/data/Method3_0067/outs/filtered_feature_bc_matrix.h5"
  ),
  fragment_paths = list(
    "0045" = "/path/to/data/Method3_0045/outs/atac_fragments.tsv.gz",
    "0051" = "/path/to/data/Method3_0051/outs/atac_fragments.tsv.gz",
    "0060" = "/path/to/data/Method3_0060/outs/atac_fragments.tsv.gz",
    "0067" = "/path/to/data/Method3_0067/outs/atac_fragments.tsv.gz"
  )
)

output_dir <- "outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# Data Import and Initial Processing
# =============================================================================

load_sample_data <- function(sample_id) {
  cat(paste("Loading sample", sample_id, "...\n"))
  
  inputdata <- Read10X_h5(sample_info$data_paths[[sample_id]])
  
  rna_counts <- inputdata$`Gene Expression`
  atac_counts <- inputdata$Peaks
  
  cat(paste("Sample", sample_id, "- RNA features:", nrow(rna_counts), 
            "| ATAC peaks:", nrow(atac_counts), "| Cells:", ncol(rna_counts), "\n"))
  
  return(list(rna = rna_counts, atac = atac_counts))
}

sample_data <- list()
for (sample_id in sample_info$samples) {
  sample_data[[sample_id]] <- load_sample_data(sample_id)
}

# =============================================================================
# Create Initial Seurat Objects with ATAC Integration
# =============================================================================

create_multiome_object <- function(sample_id, rna_counts, atac_counts) {
  cat(paste("Processing sample", sample_id, "...\n"))
  
  seurat_obj <- CreateSeuratObject(counts = rna_counts, project = paste0("M3_", sample_id))
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  atac_counts_filtered <- atac_counts[as.vector(grange_use), ]
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts_filtered,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = sample_info$fragment_paths[[sample_id]],
    min.cells = 10,
    annotation = annotations
  )
  
  seurat_obj[["ATAC"]] <- chrom_assay
  
  cat(paste("Filtered ATAC peaks for sample", sample_id, ":", nrow(atac_counts_filtered), "\n"))
  
  return(seurat_obj)
}

seurat_objects <- list()
for (sample_id in sample_info$samples) {
  seurat_objects[[sample_id]] <- create_multiome_object(
    sample_id, 
    sample_data[[sample_id]]$rna, 
    sample_data[[sample_id]]$atac
  )
}

# =============================================================================
# Quality Control Visualization
# =============================================================================

create_qc_plots <- function(seurat_obj, sample_id) {
  qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC")
  
  qc_plot <- VlnPlot(
    seurat_obj, 
    features = qc_features, 
    ncol = 5,
    log = TRUE, 
    pt.size = 0
  ) + 
    NoLegend() + 
    ggtitle(paste("Sample", sample_id, "- QC Metrics"))
  
  return(qc_plot)
}

qc_plots <- list()
for (sample_id in sample_info$samples) {
  qc_plots[[sample_id]] <- create_qc_plots(seurat_objects[[sample_id]], sample_id)
  print(qc_plots[[sample_id]])
}


# =============================================================================
# Comprehensive Analysis Pipeline
# =============================================================================

perform_complete_analysis <- function(seurat_obj, sample_id) {
  cat(paste("Analyzing sample", sample_id, "...\n"))
  
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
  
  DefaultAssay(seurat_obj) <- "ATAC"
  seurat_obj <- RunTFIDF(seurat_obj) %>%
    FindTopFeatures(min.cutoff = "q0") %>%
    RunSVD() %>%
    RunUMAP(reduction = "lsi", dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  seurat_obj <- FindMultiModalNeighbors(
    seurat_obj, 
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:20, 2:20), 
    k.nn = 10, 
    knn.range = 100
  ) %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
    FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  
  return(seurat_obj)
}

processed_objects <- list()
for (sample_id in sample_info$samples) {
  processed_objects[[sample_id]] <- perform_complete_analysis(seurat_objects[[sample_id]], sample_id)
  
  saveRDS(processed_objects[[sample_id]], 
          file.path(output_dir, paste0("M3_", sample_id, "_processed.RDS")))
}

# =============================================================================
# UMAP Visualization
# =============================================================================

create_umap_plots <- function(seurat_obj, sample_id) {
  p1 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste("RNA - Sample", sample_id))
  
  p2 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste("ATAC - Sample", sample_id))
  
  p3 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste("WNN - Sample", sample_id))
  
  combined_plot <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  
  return(combined_plot)
}


umap_plots <- list()
for (sample_id in sample_info$samples) {
  umap_plots[[sample_id]] <- create_umap_plots(processed_objects[[sample_id]], sample_id)
  print(umap_plots[[sample_id]])
}



# =============================================================================
# Save Results
# =============================================================================

for (sample_id in sample_info$samples) {
  saveRDS(processed_objects[[sample_id]], 
          file.path(output_dir, paste0("M3_", sample_id, "_WNN.RDS")))
}
