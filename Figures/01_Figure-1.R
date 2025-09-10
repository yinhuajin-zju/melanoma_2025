# =============================================================================
# Figure-1 QC 
# Author: Hanzhang
# Description: Visualization of data with focus on quality control
#              metrics across multiple samples
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scCATCH)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  library(patchwork)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(biovizBase)
  library(hdf5r)
  library(dplyr)
  library(readxl)
  library(ggsignif)
})

# =============================================================================
# Global Configuration
# =============================================================================

data_path <- "/path/to/your/data/"
output_path <- "/path/to/your/outputs/"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

custom_colors <- list(
  primary = c("#E57164", "#A184BC", "#0EACC9", "#FFB547", "#62A870"),
  extended = c("#E57164", "#A184BC", "#0EACC9", "#FFB547", "#62A870", 
               "#F48C77", "#BE9DD2", "#39BEDD", "#FFC970", "#82C28A", 
               "#CC4C3A", "#514A83")
)

marker_genes <- list(
  immune = c("CD3E", "CD3D", "CD2", "CD4", "CD8A", "CD8B"),
  nk = c("NCAM1", "KLRB1", "NKG7")
)

# =============================================================================
# Data File Configuration
# =============================================================================

data_files <- list(
  M1_fix = "M1_fix_wnn.RDS",
  M1_nofix = "M1_nofix_wnn.RDS",
  M2_fix = "M2_fix_wnn.RDS",
  M2_nofix = "M2_nofix_wnn.RDS",
  M3_integrated = "M3_integrated.RDS",
  M3_0045 = "M3_0045_wnn.RDS",
  M3_0051 = "M3_0051_wnn.RDS",
  M3_0060 = "M3_0060_wnn.RDS",
  M3_0067 = "M3_0067_wnn.RDS"
)

load_seurat_objects <- function(file_list, base_path = data_path) {
  objects <- list()
  
  for (name in names(file_list)) {
    file_path <- file.path(base_path, file_list[[name]])
    
    if (file.exists(file_path)) {
      cat(paste("Loading", name, "from", file_path, "...\n"))
      objects[[name]] <- readRDS(file_path)
      cat(paste("  - Loaded", ncol(objects[[name]]), "cells\n"))
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  
  return(objects)
}

M2_nofix_path <- file.path(data_path, data_files$M2_nofix)
if (file.exists(M2_nofix_path)) {
  M2_nofix <- readRDS(M2_nofix_path)
  cat(paste("Loaded M2_nofix with", ncol(M2_nofix), "cells\n"))
} else {
  stop("M2_nofix file not found. Please check the file path.")
}

# =============================================================================
# Quality Control Data Extraction Functions
# =============================================================================

#' Extract QC metrics from Seurat object
#' @param seurat_object Seurat object
#' @param sample_name Sample identifier (optional)
#' @return QC data frame
extract_qc_data <- function(seurat_object, sample_name = "Sample") {
  if (!inherits(seurat_object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  
  required_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  missing_cols <- setdiff(required_cols, colnames(seurat_object@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required metadata columns:", paste(missing_cols, collapse = ", ")))
  }
  
  qc_data <- data.frame(
    Cell = colnames(seurat_object),
    Sample = sample_name,
    nFeature_RNA = seurat_object$nFeature_RNA,
    nCount_RNA = seurat_object$nCount_RNA,
    percent.mt = seurat_object$percent.mt,
    stringsAsFactors = FALSE
  )
  
  if ("nFeature_ATAC" %in% colnames(seurat_object@meta.data)) {
    qc_data$nFeature_ATAC <- seurat_object$nFeature_ATAC
    qc_data$nCount_ATAC <- seurat_object$nCount_ATAC
  }
  
  return(qc_data)
}

qc_results <- extract_qc_data(M2_nofix, sample_name = "M2_nofix")
print(summary(qc_results[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]))

# =============================================================================
# Quality Control Visualization Functions
# =============================================================================

#' Create QC scatter plot
#' @param data QC data frame
#' @param title Plot title
#' @return ggplot object
plot_qc_scatter <- function(data, title = "RNA Feature vs Count Distribution") {
  if (!all(c("nFeature_RNA", "nCount_RNA", "percent.mt") %in% colnames(data))) {
    stop("Required columns not found in data.")
  }
  
  ggplot(data, aes(x = nFeature_RNA, y = nCount_RNA, color = percent.mt)) +
    geom_point(alpha = 0.6, size = 3) +
    scale_color_gradient(
      low = custom_colors$primary[3], 
      high = custom_colors$primary[1],
      name = "Percent\nMitochondrial"
    ) +
    labs(
      title = title,
      subtitle = "Colored by mitochondrial gene percentage",
      x = "Number of RNA Features", 
      y = "RNA Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )
}

#' Create grouped points plot by mitochondrial percentage
#' @param data QC data frame
#' @param title Plot title
#' @return ggplot object
plot_qc_grouped_points <- function(data, title = "RNA Features by Mitochondrial Percentage Groups") {
  # Check required columns
  if (!all(c("percent.mt", "nFeature_RNA") %in% colnames(data))) {
    stop("Required columns 'percent.mt' and 'nFeature_RNA' not found.")
  }
  
  data$percent.mt.group <- cut(
    data$percent.mt, 
    breaks = c(0, 10, 20, 30, 50, 100),
    labels = c("0-10%", "10-20%", "20-30%", "30-50%", ">50%"),
    include.lowest = TRUE
  )
  
  ggplot(data, aes(x = percent.mt.group, y = nFeature_RNA, color = percent.mt.group)) +
    geom_point(
      alpha = 0.6, 
      position = position_jitter(width = 0.2), 
      size = 2.5, 
      stroke = 1.2
    ) +
    stat_summary(
      fun = median, 
      geom = "crossbar", 
      width = 0.5, 
      color = "black", 
      size = 0.8
    ) +
    scale_color_manual(
      values = custom_colors$primary, 
      name = "Mitochondrial\nPercentage"
    ) +
    labs(
      title = title,
      subtitle = "Points show individual cells, crossbars show medians",
      x = "Mitochondrial Percentage Group", 
      y = "Number of RNA Features"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )
}

#' Create pie chart for mitochondrial percentage distribution
#' @param data QC data frame
#' @param title Plot title
#' @return ggplot object
plot_mt_distribution_pie <- function(data, title = "Distribution of Mitochondrial Gene Percentage") {

    data$percent.mt.group <- cut(
    data$percent.mt, 
    breaks = c(0, 10, 20, 30, 50, 100),
    labels = c("0-10%", "10-20%", "20-30%", "30-50%", ">50%"),
    include.lowest = TRUE
  )
  
  # Calculate percentages
  pie_data <- data %>%
    group_by(percent.mt.group) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(
      percentage = round(count / sum(count) * 100, 1),
      label = paste0(percentage, "%")
    )
  
  ggplot(pie_data, aes(x = "", y = percentage, fill = percent.mt.group)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
    coord_polar("y", start = 0) +
    geom_text(
      aes(label = label), 
      position = position_stack(vjust = 0.5), 
      size = 5, 
      fontface = "bold", 
      color = "white"
    ) +
    scale_fill_manual(
      values = custom_colors$primary, 
      name = "Mitochondrial\nPercentage"
    ) +
    labs(
      title = title,
      subtitle = "Cell proportion by mitochondrial content"
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
}

# =============================================================================
# Generate Quality Control Visualizations
# =============================================================================

qc_scatter <- plot_qc_scatter(qc_results, "M2_nofix: RNA Feature vs Count Distribution")
qc_grouped <- plot_qc_grouped_points(qc_results, "M2_nofix: RNA Features by Mitochondrial Groups")
qc_pie <- plot_mt_distribution_pie(qc_results, "M2_nofix: Mitochondrial Gene Distribution")

# Display plots
print(qc_scatter)
print(qc_grouped)
print(qc_pie)

# =============================================================================
# Save Visualization Results
# =============================================================================

ggsave(
  filename = file.path(output_path, "qc_scatter_plot.pdf"),
  plot = qc_scatter,
  width = 10, height = 8
)

ggsave(
  filename = file.path(output_path, "qc_grouped_points.pdf"),
  plot = qc_grouped,
  width = 12, height = 8
)

ggsave(
  filename = file.path(output_path, "qc_pie_chart.pdf"),
  plot = qc_pie,
  width = 10, height = 8
)

combined_qc_plot <- (qc_scatter | qc_pie) / qc_grouped +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Comprehensive Quality Control Analysis",
    subtitle = "Multi-panel view of RNA and mitochondrial gene metrics",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

print(combined_qc_plot)

ggsave(
  filename = file.path(output_path, "combined_qc_analysis.pdf"),
  plot = combined_qc_plot,
  width = 20, height = 15
)
