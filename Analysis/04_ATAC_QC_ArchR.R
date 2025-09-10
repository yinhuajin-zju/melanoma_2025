# =============================================================================
# ArchR analysis Method 2: Fixed vs Non-fixed Comparison
# Author: Hanzhang
# Description: Comprehensive ATAC-seq analysis comparing Method 2 fixed and 
#              non-fixed samples using ArchR for quality control, doublet 
#              detection, and visualization
# =============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(pheatmap)
  library(Rsamtools)
  library(scran)
  library(scater)
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(SingleCellExperiment)
  library(ComplexHeatmap)
  library(ggplot2)
  library(stringr)
  library(EnsDb.Hsapiens.v86)
  library(viridis)
})

# =============================================================================
# ArchR Configuration and Setup
# =============================================================================

addArchRGenome("hg38")  # Human genome hg38
# addArchRGenome("mm10") # Uncomment for mouse genome
addArchRThreads(threads = 10)  # Adjust based on available cores

# =============================================================================
# File Path Configuration
# =============================================================================

setwd("/path/to/your/method2/directory")

input_fragment_files <- c(
  "/path/to/data/Method2_nofix/outs/atac_fragments.tsv.gz",
  "/path/to/data/Method2_fix/outs/atac_fragments.tsv.gz"
)

sample_names <- c(
  "Method2_nofix",
  "Method2_fix"
)

for (i in seq_along(input_fragment_files)) {
  if (file.exists(input_fragment_files[i])) {
    cat(paste("✓ Found:", sample_names[i], "\n"))
  } else {
    stop(paste("✗ File not found:", input_fragment_files[i]))
  }
}

output_dir <- "Method2_ATAC_analysis_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# Arrow File Creation
# =============================================================================

ArrowFiles <- createArrowFiles(
  inputFiles = input_fragment_files,
  sampleNames = sample_names,
  minTSS = 4,        # Minimum TSS enrichment score
  minFrags = 1000,   # Minimum number of fragments per cell
  addTileMat = TRUE, # Add tile matrix for peak calling
  addGeneScoreMat = TRUE, # Add gene score matrix
  excludeChr = c("chrM", "chrY", "chrX") # Exclude problematic chromosomes
)

# =============================================================================
# Doublet Detection
# =============================================================================

doublet_scores <- addDoubletScores(
  input = ArrowFiles,
  k = 20,           # Number of cells near pseudo-doublet to count
  knnMethod = "UMAP",
  useMatrix = "TileMatrix",
  nTrials = 5,      # Number of trials for doublet simulation
  LSIMethod = 1,
  scaleDims = FALSE,
  corCutOff = 0.75, # Correlation cutoff for doublet identification
  UMAPParams = list(
    n_neighbors = 30,
    min_dist = 0.3,
    metric = "cosine",
    verbose = TRUE
  )
)

# =============================================================================
# ArchR Project Creation
# =============================================================================

proj_initial <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = output_dir,
  copyArrows = TRUE  # Create copies of Arrow files for safety
)

cat("Initial cell counts by sample:\n")
initial_counts <- table(proj_initial@cellColData$Sample)
print(initial_counts)

# =============================================================================
# Quality Control Filtering
# =============================================================================

proj_filtered <- filterDoublets(proj_initial)

cat("Cell counts after doublet filtering:\n")
filtered_counts <- table(proj_filtered@cellColData$Sample)
print(filtered_counts)

total_initial <- sum(initial_counts)
total_filtered <- sum(filtered_counts)
percent_retained <- round((total_filtered / total_initial) * 100, 2)

cat(paste("Filtering summary:\n"))
cat(paste("  Initial cells:", total_initial, "\n"))
cat(paste("  Filtered cells:", total_filtered, "\n"))
cat(paste("  Retention rate:", percent_retained, "%\n"))

# =============================================================================
# Quality Control Visualization
# =============================================================================

fragment_plot <- plotFragmentSizes(ArchRProj = proj_filtered) +
  ggtitle("Method 2: Fragment Size Distribution") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

tss_plot <- plotTSSEnrichment(ArchRProj = proj_filtered) +
  ggtitle("Method 2: TSS Enrichment Profile") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

plotPDF(
  fragment_plot, 
  tss_plot, 
  name = "Method2-QC-FragSizes-TSSProfile_filtered.pdf", 
  ArchRProj = proj_filtered, 
  addDOC = FALSE, 
  width = 5, 
  height = 5
)

# =============================================================================
# Sample-specific Density Plots
# =============================================================================

filtered_sample_counts <- table(proj_filtered$Sample)

sample_density_plots <- list()

for (i in seq_along(sample_names)) {
  cat(paste("Processing", sample_names[i], "...\n"))
  
  proj_sample <- proj_filtered[proj_filtered$Sample == sample_names[i]]
  
  density_plot <- ggPoint(
    x = log10(proj_sample$nFrags),
    y = proj_sample$TSSEnrichment,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10(Unique Fragments)",
    ylabel = "TSS Enrichment"
  ) + 
    geom_hline(yintercept = 4, lty = "dashed", color = "red", size = 1) + 
    geom_vline(xintercept = 3, lty = "dashed", color = "red", size = 1) +
    ggtitle(paste0(sample_names[i], "\n", "Cells Pass Filter = ", 
                   filtered_sample_counts[[sample_names[i]]])) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  sample_density_plots[[i]] <- density_plot
}

plotPDF(
  sample_density_plots[[1]], 
  sample_density_plots[[2]],
  name = "Method2-QC-TSSenrich-Fragment_density_comparison.pdf", 
  ArchRProj = proj_filtered, 
  addDOC = FALSE, 
  width = 5, 
  height = 5
)

# =============================================================================
# Fixed vs Non-fixed Comparison Visualization
# =============================================================================

comparison_plot <- sample_density_plots[[1]] + sample_density_plots[[2]] +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Method 2: Fixed vs Non-fixed ATAC-seq Quality Comparison",
    subtitle = "TSS Enrichment vs Fragment Count Density",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  filename = file.path(output_dir, "Method2_fixed_vs_nofix_comparison.pdf"),
  plot = comparison_plot,
  width = 12, height = 6
)

# =============================================================================
# Save Project and Generate Summary
# =============================================================================

saveArchRProject(ArchRProj = proj_filtered, outputDirectory = output_dir, load = FALSE)

