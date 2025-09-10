# =============================================================================
# Method 2 Analysis: Comparison of Fixed vs Non-fixed Samples
# Author: Hanzhang
# Description: Data load and analysis of scRNA-seq and scATAC-seq data from 
#              fixed and non-fixed Method 2 samples using Seurat and Signac 
#              workflows
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
# Configuration and File Paths
# =============================================================================

setwd("/path/to/your/method2/directory")

data_paths <- list(
  nofix = "/path/to/data/Method2_nofix/outs/filtered_feature_bc_matrix.h5",
  fix = "/path/to/data/Method2_fix/outs/filtered_feature_bc_matrix.h5"
)

fragment_paths <- list(
  nofix = "/path/to/data/Method2_nofix/outs/atac_fragments.tsv.gz",
  fix = "/path/to/data/Method2_fix/outs/atac_fragments.tsv.gz"
)

output_dir <- "outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# Data Import and Initial Processing
# =============================================================================

inputdata.10x_nofix <- Read10X_h5(data_paths$nofix)
inputdata.10x_fix <- Read10X_h5(data_paths$fix)

rna_counts_nofix <- inputdata.10x_nofix$`Gene Expression`
atac_counts_nofix <- inputdata.10x_nofix$Peaks

rna_counts_fix <- inputdata.10x_fix$`Gene Expression`
atac_counts_fix <- inputdata.10x_fix$Peaks

cat(paste("RNA features (nofix):", nrow(rna_counts_nofix), "| Cells:", ncol(rna_counts_nofix), "\n"))
cat(paste("RNA features (fix):", nrow(rna_counts_fix), "| Cells:", ncol(rna_counts_fix), "\n"))
cat(paste("ATAC peaks (nofix):", nrow(atac_counts_nofix), "| Cells:", ncol(atac_counts_nofix), "\n"))
cat(paste("ATAC peaks (fix):", nrow(atac_counts_fix), "| Cells:", ncol(atac_counts_fix), "\n"))

# =============================================================================
# Create Initial Seurat Objects
# =============================================================================

M2_fix <- CreateSeuratObject(counts = rna_counts_fix, project = "M2_fix")
M2_nofix <- CreateSeuratObject(counts = rna_counts_nofix, project = "M2_nofix")

M2_fix[["percent.mt"]] <- PercentageFeatureSet(M2_fix, pattern = "^MT-")
M2_nofix[["percent.mt"]] <- PercentageFeatureSet(M2_nofix, pattern = "^MT-")

cat("M2_fix initial metadata:\n")
print(head(M2_fix@meta.data))

cat("M2_nofix initial metadata:\n")
print(head(M2_nofix@meta.data))

# =============================================================================
# Initial RNA Quality Control Visualization
# =============================================================================

rna_qc_features <- c("percent.mt", "nCount_RNA", "nFeature_RNA")

qc_plot_fix <- VlnPlot(
  M2_fix, 
  features = rna_qc_features, 
  pt.size = 0, 
  cols = "red"
) + 
  ggtitle("M2_fix - Initial RNA QC")

qc_plot_nofix <- VlnPlot(
  M2_nofix, 
  features = rna_qc_features, 
  pt.size = 0, 
  cols = "red"
) + 
  ggtitle("M2_nofix - Initial RNA QC")

print(qc_plot_fix)
print(qc_plot_nofix)

# =============================================================================
# ATAC Data Processing
# =============================================================================

process_atac_data <- function(atac_counts, sample_name) {
  cat(paste("Processing ATAC data for", sample_name, "...\n"))
  
  grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  atac_counts_filtered <- atac_counts[as.vector(grange_use), ]
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  
  cat(paste("Filtered ATAC peaks for", sample_name, ":", nrow(atac_counts_filtered), "\n"))
  
  return(list(counts = atac_counts_filtered, annotations = annotations))
}

atac_processed_nofix <- process_atac_data(atac_counts_nofix, "M2_nofix")
atac_processed_fix <- process_atac_data(atac_counts_fix, "M2_fix")

# =============================================================================
# Create and Add Chromatin Assays
# =============================================================================

create_chromatin_assay <- function(atac_counts, fragment_file, annotations, sample_name) {
  cat(paste("Creating chromatin assay for", sample_name, "...\n"))
  
  chromatin_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = fragment_file,
    min.cells = 10,
    annotation = annotations
  )
  
  return(chromatin_assay)
}

M2_fix[["ATAC"]] <- create_chromatin_assay(
  atac_processed_fix$counts,
  fragment_paths$fix,
  atac_processed_fix$annotations,
  "M2_fix"
)

M2_nofix[["ATAC"]] <- create_chromatin_assay(
  atac_processed_nofix$counts,
  fragment_paths$nofix,
  atac_processed_nofix$annotations,
  "M2_nofix"
)

# =============================================================================
# Comprehensive Quality Control Visualization
# =============================================================================

all_qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC")

comprehensive_qc_fix <- VlnPlot(
  M2_fix,
  features = all_qc_features,
  ncol = 5,
  log = TRUE,
  pt.size = 0
) + 
  NoLegend() + 
  ggtitle("M2_fix - Comprehensive QC Metrics")

comprehensive_qc_nofix <- VlnPlot(
  M2_nofix,
  features = all_qc_features,
  ncol = 5,
  log = TRUE,
  pt.size = 0
) + 
  NoLegend() + 
  ggtitle("M2_nofix - Comprehensive QC Metrics")

print(comprehensive_qc_fix)
print(comprehensive_qc_nofix)

# =============================================================================
# Quality Control Filtering
# =============================================================================

filtering_criteria <- list(
  nCount_ATAC = c(5e3, 7e4),
  nCount_RNA = c(1000, 22000),
  percent.mt = c(0, 20)
)

cat(paste("Pre-filtering - M2_fix cells:", ncol(M2_fix), "\n"))
cat(paste("Pre-filtering - M2_nofix cells:", ncol(M2_nofix), "\n"))

M2_fix <- subset(
  x = M2_fix,
  subset = nCount_ATAC > filtering_criteria$nCount_ATAC[1] &
    nCount_ATAC < filtering_criteria$nCount_ATAC[2] &
    nCount_RNA > filtering_criteria$nCount_RNA[1] &
    nCount_RNA < filtering_criteria$nCount_RNA[2] &
    percent.mt < filtering_criteria$percent.mt[2]
)

M2_nofix <- subset(
  x = M2_nofix,
  subset = nCount_ATAC > filtering_criteria$nCount_ATAC[1] &
    nCount_ATAC < filtering_criteria$nCount_ATAC[2] &
    nCount_RNA > filtering_criteria$nCount_RNA[1] &
    nCount_RNA < filtering_criteria$nCount_RNA[2] &
    percent.mt < filtering_criteria$percent.mt[2]
)

cat(paste("Post-filtering - M2_fix cells:", ncol(M2_fix), "\n"))
cat(paste("Post-filtering - M2_nofix cells:", ncol(M2_nofix), "\n"))

# =============================================================================
# RNA Data Processing and Dimension Reduction
# =============================================================================

process_rna_workflow <- function(seurat_obj, sample_name) {
  cat(paste("Processing RNA data for", sample_name, "...\n"))
  
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
  
  return(seurat_obj)
}

M2_fix <- process_rna_workflow(M2_fix, "M2_fix")
M2_nofix <- process_rna_workflow(M2_nofix, "M2_nofix")

# =============================================================================
# ATAC Data Processing and Dimension Reduction
# =============================================================================

process_atac_workflow <- function(seurat_obj, sample_name) {
  cat(paste("Processing ATAC data for", sample_name, "...\n"))
  
  DefaultAssay(seurat_obj) <- "ATAC"
  seurat_obj <- RunTFIDF(seurat_obj) %>%
    FindTopFeatures(min.cutoff = "q0") %>%
    RunSVD() %>%
    RunUMAP(reduction = "lsi", dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  return(seurat_obj)
}

M2_fix <- process_atac_workflow(M2_fix, "M2_fix")
M2_nofix <- process_atac_workflow(M2_nofix, "M2_nofix")

# =============================================================================
# Weighted Nearest Neighbor (WNN) Analysis
# =============================================================================

perform_wnn_workflow <- function(seurat_obj, sample_name, use_custom_params = FALSE) {
  cat(paste("Running WNN analysis for", sample_name, "...\n"))
  
  if (use_custom_params) {
    seurat_obj <- FindMultiModalNeighbors(
      seurat_obj,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:20, 2:20),
      k.nn = 10,
      knn.range = 100
    )
  } else {
    seurat_obj <- FindMultiModalNeighbors(
      seurat_obj,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:20, 2:20)
    )
  }
  
  seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
    FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  
  return(seurat_obj)
}

M2_nofix <- perform_wnn_workflow(M2_nofix, "M2_nofix", use_custom_params = FALSE)
M2_fix <- perform_wnn_workflow(M2_fix, "M2_fix", use_custom_params = TRUE)

# =============================================================================
# Visualization: UMAP Plots
# =============================================================================

create_umap_plots <- function(seurat_obj, sample_name) {
  p1 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste(sample_name, "- RNA"))
  
  p2 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste(sample_name, "- ATAC"))
  
  p3 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "seurat_clusters", 
                label = TRUE, label.size = 2.5, repel = TRUE) + 
    ggtitle(paste(sample_name, "- WNN"))
  
  combined_plot <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  
  return(combined_plot)
}

umap_plots_nofix <- create_umap_plots(M2_nofix, "M2_nofix")
umap_plots_fix <- create_umap_plots(M2_fix, "M2_fix")

print(umap_plots_nofix)
print(umap_plots_fix)

# =============================================================================
# Differential Expression Analysis
# =============================================================================

find_all_markers_safe <- function(seurat_obj, sample_name) {
  cat(paste("Finding markers for", sample_name, "...\n"))
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj[["RNA"]]$data <- as(object = seurat_obj[["RNA"]]$counts, Class = "dgCMatrix")
  
  all_markers <- FindAllMarkers(object = seurat_obj, verbose = TRUE)
  
  return(all_markers)
}

all_markers_M2_nofix <- find_all_markers_safe(M2_nofix, "M2_nofix")
all_markers_M2_fix <- find_all_markers_safe(M2_fix, "M2_fix")

cat("Top markers for M2_nofix:\n")
print(head(all_markers_M2_nofix))

cat("Top markers for M2_fix:\n")
print(head(all_markers_M2_fix))

# =============================================================================
# Feature Visualization
# =============================================================================

genes_of_interest <- c("NCAM1", "PTPRC")

# Function to create feature plots
create_feature_plots <- function(seurat_obj, genes, sample_name) {
  plot_list <- list()
  
  for (gene in genes) {
    plot_list[[paste(sample_name, gene, sep = "_")]] <- FeaturePlot(
      object = seurat_obj,
      features = gene,
      reduction = "wnn.umap",
      min.cutoff = "q9",
      cols = c("lightgrey", "blue")
    ) + 
      ggtitle(paste(sample_name, "-", gene))
  }
  
  return(plot_list)
}

feature_plots_fix <- create_feature_plots(M2_fix, genes_of_interest, "M2_fix")
feature_plots_nofix <- create_feature_plots(M2_nofix, genes_of_interest, "M2_nofix")

for (plot_name in names(feature_plots_fix)) {
  print(feature_plots_fix[[plot_name]])
}

for (plot_name in names(feature_plots_nofix)) {
  print(feature_plots_nofix[[plot_name]])
}

cat("Creating comparison feature plots...\n")

ncam1_comparison <- feature_plots_fix[["M2_fix_NCAM1"]] + feature_plots_nofix[["M2_nofix_NCAM1"]] +
  plot_annotation(title = "NCAM1 Expression Comparison")

ptprc_comparison <- feature_plots_fix[["M2_fix_PTPRC"]] + feature_plots_nofix[["M2_nofix_PTPRC"]] +
  plot_annotation(title = "PTPRC Expression Comparison")

print(ncam1_comparison)
print(ptprc_comparison)

# =============================================================================
# Save Results
# =============================================================================

saveRDS(M2_fix, file.path(output_dir, "M2_fix_wnn.RDS"))
saveRDS(M2_nofix, file.path(output_dir, "M2_nofix_wnn.RDS"))


