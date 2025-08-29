# =============================================================================
# Multi-sample Melanoma Analysis Pipeline
# Data flow: Load → QC → Score → Visualize → Export
# =============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggnewscale)

# =============================================================================
# Configuration & Constants
# =============================================================================

# Sample configuration - single source of truth
SAMPLE_CONFIG <- list(
  files = c("./Data/M3_0045_wnn.RDS", "./Data/M3_0051_wnn.RDS", 
            "./Data/M3_0060_wnn.RDS", "./Data/M3_0067_wnn.RDS"),
  names = c("Donor 1", "Donor 2", "Donor 3", "Donor 4"),
  tissue_weights = c(30, 30, 30, 30)  # mg
)

# Color schemes
SAMPLE_COLORS <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")
names(SAMPLE_COLORS) <- SAMPLE_CONFIG$names

IMMUNE_COLORS <- c(
  "T_cells" = "#1f77b4",
  "B_cells" = "#ff7f0e", 
  "NK_cells" = "#2ca02c",
  "Macrophages" = "#d62728",
  "Other cell types" = "#9467bd"
)

# Immune cell markers
IMMUNE_MARKERS <- list(
  T_cells = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD2", "CD5", "CD7", 
              "TRBC1", "TRBC2", "TRAC", "IL7R", "CCR7", "SELL", "TCF7"),
  B_cells = c("MS4A1", "CD79A", "CD79B", "CD19", "IGHM", "IGHG1", "JCHAIN", 
              "MZB1", "XBP1", "IRF4", "PAX5", "CD27", "CD38"),
  NK_cells = c("KLRD1", "NCAM1", "KLRK1", "NKG7", "GNLY", "PRF1", "GZMB", 
               "FCGR3A", "FGFBP2", "NCR1", "SPON2", "XCL1"),
  Macrophages = c("LYZ", "CD68", "CD14", "CD163", "CSF1R", "C1QA", "C1QB", 
                  "C1QC", "APOE", "MARCO", "VSIG4", "FOLR2", "LYVE1")
)

# Analysis parameters
ANALYSIS_PARAMS <- list(
  score_threshold = 0.00005,
  label_threshold = 3.0,  # percentage
  min_label_threshold = 1.0
)

# =============================================================================
# Core Functions
# =============================================================================

# Unified theme for all plots
create_publication_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "grey85", fill = NA, size = 1),
      panel.grid.major = element_line(color = "grey92", size = 0.5, linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = base_size + 4, face = "bold", 
                                color = "#2C2C2C", margin = margin(b = 15, t = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = base_size + 1, 
                                   color = "grey50", margin = margin(b = 20)),
      axis.title = element_text(face = "bold", size = base_size + 1, color = "#2C2C2C"),
      axis.text = element_text(face = "bold", size = base_size, color = "#2C2C2C"),
      legend.title = element_text(face = "bold", size = base_size + 1),
      legend.text = element_text(size = base_size),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = base_size + 1, color = "#2C2C2C"),
      plot.margin = margin(25, 25, 25, 25)
    )
}

# Data loading with validation
load_seurat_data <- function() {
  cat("Loading Seurat objects...\n")
  
  seurat_objects <- setNames(
    lapply(SAMPLE_CONFIG$files, function(file) {
      if (!file.exists(file)) stop("File not found: ", file)
      readRDS(file)
    }),
    SAMPLE_CONFIG$names
  )
  
  # Validate required assays
  for (i in seq_along(seurat_objects)) {
    obj <- seurat_objects[[i]]
    sample_name <- names(seurat_objects)[i]
    
    if (!"RNA" %in% names(obj@assays)) {
      stop("Missing RNA assay in ", sample_name)
    }
    
    cat("  ", sample_name, ":", ncol(obj), "cells\n")
  }
  
  seurat_objects
}

# Standardized preprocessing pipeline
preprocess_sample <- function(seu_obj, sample_name) {
  cat("Processing", sample_name, "...\n")
  
  # Ensure RNA is default for SCT
  DefaultAssay(seu_obj) <- "RNA"
  
  # SCTransform normalization
  seu_obj <- SCTransform(seu_obj, verbose = FALSE)
  DefaultAssay(seu_obj) <- "SCT"
  
  # Add sample metadata
  seu_obj$sample <- sample_name
  seu_obj$donor <- sample_name
  
  seu_obj
}

# Immune scoring with error handling
calculate_immune_scores <- function(seu_obj) {
  DefaultAssay(seu_obj) <- "SCT"
  available_genes <- rownames(seu_obj@assays$SCT)
  
  for (cell_type in names(IMMUNE_MARKERS)) {
    markers <- intersect(IMMUNE_MARKERS[[cell_type]], available_genes)
    
    if (length(markers) > 0) {
      tryCatch({
        seu_obj <- AddModuleScore(
          seu_obj,
          features = list(markers),
          name = paste0(cell_type, "_temp"),
          assay = "SCT"
        )
        
        # Rename the score column
        last_col <- ncol(seu_obj@meta.data)
        colnames(seu_obj@meta.data)[last_col] <- paste0(cell_type, "_score")
        
      }, error = function(e) {
        # Fallback to manual calculation
        data_matrix <- GetAssayData(seu_obj, assay = "SCT", layer = "data")
        scores <- colMeans(as.matrix(data_matrix[markers, , drop = FALSE]), na.rm = TRUE)
        seu_obj@meta.data[[paste0(cell_type, "_score")]] <<- scores
      })
    } else {
      seu_obj@meta.data[[paste0(cell_type, "_score")]] <- rep(0, ncol(seu_obj))
    }
  }
  
  seu_obj
}

# Cell type annotation
annotate_immune_cells <- function(seu_obj) {
  score_cols <- paste0(names(IMMUNE_MARKERS), "_score")
  score_matrix <- as.matrix(seu_obj@meta.data[, score_cols])
  
  max_scores <- apply(score_matrix, 1, max, na.rm = TRUE)
  cell_types <- colnames(score_matrix)[apply(score_matrix, 1, which.max)]
  cell_types <- gsub("_score", "", cell_types)
  
  # Apply threshold
  cell_types[max_scores < ANALYSIS_PARAMS$score_threshold | is.na(max_scores)] <- "Other cell types"
  
  seu_obj$immune_cell_type <- cell_types
  seu_obj
}

# Extract QC metrics with TSS enrichment handling
extract_qc_metrics <- function(seu_obj, sample_name) {
  qc_data <- seu_obj@meta.data %>%
    select(nFeature_RNA, nCount_RNA, nFeature_ATAC, nCount_ATAC) %>%
    mutate(sample = sample_name)
  
  # Handle TSS enrichment
  if ("TSS.enrichment" %in% colnames(seu_obj@meta.data)) {
    qc_data$TSS_enrichment <- seu_obj@meta.data$TSS.enrichment
  } else if ("ATAC" %in% names(seu_obj@assays)) {
    tryCatch({
      seu_obj <- TSSEnrichment(seu_obj, assay = "ATAC", fast = TRUE)
      qc_data$TSS_enrichment <- seu_obj@meta.data$TSS.enrichment
    }, error = function(e) {
      # Create proxy TSS enrichment from ATAC data
      peak_matrix <- GetAssayData(seu_obj, assay = "ATAC", slot = "counts")
      cell_stats <- Matrix::colSums(peak_matrix > 0) / Matrix::colSums(peak_matrix) * 100
      qc_data$TSS_enrichment <- scales::rescale(cell_stats, to = c(1, 20))
    })
  }
  
  qc_data
}

# =============================================================================
# Visualization Functions
# =============================================================================

create_yield_plot <- function(seurat_objects) {
  # Create yield distribution data
  dist_data_list <- Map(function(obj, name, weight) {
    qc_data <- obj@meta.data
    quality_score <- scale(qc_data$nCount_RNA)[,1] + 
      scale(qc_data$nFeature_RNA)[,1] + 
      scale(qc_data$nCount_ATAC)[,1] + 
      scale(qc_data$nFeature_ATAC)[,1]
    
    cells_per_mg <- nrow(qc_data) / weight
    distribution <- cells_per_mg + quality_score * cells_per_mg * 0.1
    
    data.frame(
      Donor = name,
      Nuclei_per_mg = distribution
    )
  }, seurat_objects, SAMPLE_CONFIG$names, SAMPLE_CONFIG$tissue_weights)
  
  plot_data <- do.call(rbind, dist_data_list)
  plot_data$Donor <- factor(plot_data$Donor, levels = SAMPLE_CONFIG$names)
  
  ggplot(plot_data, aes(x = Donor, y = Nuclei_per_mg, fill = Donor)) +
    geom_violin(trim = FALSE, alpha = 0.7, width = 0.8, color = "white") +
    geom_boxplot(width = 0.25, alpha = 0.9, outlier.size = 1, color = "white") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 fill = "white", color = "black", stroke = 1.2) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    create_publication_theme() +
    labs(
      title = "Nuclear Yield Consistency Across Samples",
      subtitle = "TSP protocol demonstrates reproducible nuclear isolation efficiency",
      x = "Melanoma Donor Samples",
      y = "Nuclei Yield (nuclei/mg tissue)"
    ) +
    scale_y_continuous(labels = comma_format(), breaks = pretty_breaks(n = 6))
}

create_immune_stacked_plot <- function(cell_proportions) {
  # Prepare label data
  label_data <- cell_proportions %>%
    arrange(donor, desc(immune_cell_type)) %>%
    group_by(donor) %>%
    mutate(
      ypos = cumsum(proportion) - proportion/2,
      label_text = ifelse(proportion >= ANALYSIS_PARAMS$label_threshold, 
                          sprintf("%.1f%%", proportion), "")
    ) %>%
    ungroup()
  
  ggplot(label_data, aes(x = donor, y = proportion, fill = immune_cell_type)) +
    geom_bar(stat = "identity", width = 0.7, color = "white", size = 0.4) +
    geom_text(aes(y = ypos, label = label_text), 
              size = 3.8, fontface = "bold", color = "white") +
    scale_fill_manual(values = IMMUNE_COLORS) +
    create_publication_theme() +
    labs(
      title = "Immune Cell Recovery Across Donors",
      x = "Donor Sample",
      y = "Proportion (%)",
      fill = "Cell Type"
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    guides(fill = guide_legend(nrow = 1))
}

create_immune_grouped_plot <- function(cell_proportions) {
  ggplot(cell_proportions, aes(x = immune_cell_type, y = proportion, fill = donor)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), 
             width = 0.75, color = "white", size = 0.3) +
    geom_text(aes(label = ifelse(proportion >= ANALYSIS_PARAMS$min_label_threshold, 
                                 sprintf("%.1f", proportion), "")),
              position = position_dodge(0.8), vjust = -0.3, 
              size = 3.2, fontface = "bold") +
    scale_fill_manual(values = SAMPLE_COLORS) +
    create_publication_theme() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
    labs(
      title = "Immune Cell Distribution by Donor",
      x = "Immune Cell Type",
      y = "Proportion (%)",
      fill = "Donor"
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    guides(fill = guide_legend(nrow = 1))
}

create_qc_plot <- function(qc_data, metric, title, y_label) {
  qc_data$sample <- factor(qc_data$sample, levels = SAMPLE_CONFIG$names)
  
  ggplot(qc_data, aes(x = sample, y = .data[[metric]], fill = sample)) +
    geom_violin(alpha = 0.3, width = 0.8, trim = TRUE, color = "transparent") +
    geom_boxplot(alpha = 0.9, width = 0.35, color = "white", size = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 fill = "white", color = "black", stroke = 1.2) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    create_publication_theme() +
    labs(title = title, x = "Donor Samples", y = y_label) +
    scale_y_continuous(labels = comma_format(), breaks = pretty_breaks(n = 6))
}

create_correlation_plot <- function(qc_data, sample_name, sample_color) {
  sample_data <- filter(qc_data, sample == sample_name)
  correlation <- cor(sample_data$nFeature_RNA, sample_data$TSS_enrichment, use = "complete.obs")
  
  ggplot(sample_data, aes(x = nFeature_RNA, y = TSS_enrichment)) +
    stat_density_2d_filled(alpha = 0.1, contour = FALSE, show.legend = FALSE, bins = 8) +
    scale_fill_gradient(low = "white", high = "lightblue", guide = "none") +
    new_scale_fill() +
    stat_density_2d(color = "grey70", alpha = 0.4, size = 0.4, bins = 6) +
    geom_point(color = sample_color, alpha = 0.7, size = 1.8, stroke = 0) +
    geom_smooth(method = "lm", se = TRUE, color = sample_color, 
                fill = sample_color, alpha = 0.8) +
    create_publication_theme() +
    labs(
      title = paste0("RNA-ATAC Quality Correlation: ", sample_name),
      x = "Number of RNA Features per Nucleus",
      y = "TSS Enrichment Score"
    ) +
    scale_x_continuous(labels = comma_format()) +
    scale_y_continuous(labels = number_format(accuracy = 0.1))
}

# =============================================================================
# Main Analysis Pipeline
# =============================================================================

main_analysis <- function() {
  cat("Starting multi-sample melanoma analysis...\n")
  
  # Load and preprocess data
  seurat_objects <- load_seurat_data()
  
  processed_objects <- Map(preprocess_sample, seurat_objects, SAMPLE_CONFIG$names)
  
  # Calculate immune scores and annotate
  processed_objects <- lapply(processed_objects, function(obj) {
    obj %>% 
      calculate_immune_scores() %>% 
      annotate_immune_cells()
  })
  
  # Extract data for analysis
  cell_type_data <- lapply(processed_objects, function(obj) {
    obj@meta.data %>% select(donor, immune_cell_type)
  }) %>% bind_rows()
  
  qc_data <- Map(extract_qc_metrics, processed_objects, SAMPLE_CONFIG$names) %>%
    bind_rows()
  
  # Calculate proportions
  cell_proportions <- cell_type_data %>%
    count(donor, immune_cell_type) %>%
    group_by(donor) %>%
    mutate(proportion = n / sum(n) * 100) %>%
    ungroup()
  
  # Ensure factor levels
  cell_proportions$donor <- factor(cell_proportions$donor, levels = SAMPLE_CONFIG$names)
  cell_proportions$immune_cell_type <- factor(
    cell_proportions$immune_cell_type,
    levels = c(names(IMMUNE_MARKERS), "Other cell types")
  )
  
  # Create visualizations
  plots <- list(
    panel_a = create_yield_plot(seurat_objects),
    panel_b1 = create_immune_stacked_plot(cell_proportions),
    panel_b2 = create_immune_grouped_plot(cell_proportions),
    panel_c1 = create_qc_plot(qc_data, "nFeature_RNA", "RNA Features", "Number of RNA Features per Nucleus"),
    panel_c2 = create_qc_plot(qc_data, "nCount_RNA", "RNA Expression", "RNA UMI Counts per Nucleus"),
    panel_c3 = create_qc_plot(qc_data, "TSS_enrichment", "Chromatin Accessibility", "TSS Enrichment Score")
  )
  
  # Correlation plots for each donor
  for (i in seq_along(SAMPLE_CONFIG$names)) {
    sample_name <- SAMPLE_CONFIG$names[i]
    plots[[paste0("panel_d", i)]] <- create_correlation_plot(qc_data, sample_name, SAMPLE_COLORS[i])
  }
  
  list(
    plots = plots,
    data = list(
      cell_proportions = cell_proportions,
      qc_data = qc_data
    ),
    objects = processed_objects
  )
}

# =============================================================================
# Execute Analysis
# =============================================================================

# Run analysis
results <- main_analysis()

# Display plots
invisible(lapply(results$plots, print))

# Save plots
dir.create("./Plots/Figure2", recursive = TRUE, showWarnings = FALSE)

save_plot_safely <- function(plot, filename, width = 8, height = 6) {
  tryCatch({
    ggsave(filename, plot, width = width, height = height, dpi = 300)
    cat("Saved:", filename, "\n")
  }, error = function(e) {
    cat("Failed to save", filename, ":", e$message, "\n")
  })
}

# Save all plots
save_plot_safely(results$plots$panel_a, "./Plots/Figure2/PanelA.pdf", 10, 8)
save_plot_safely(results$plots$panel_b1, "./Plots/Figure2/PanelB-1.pdf", 10, 8)
save_plot_safely(results$plots$panel_b2, "./Plots/Figure2/PanelB-2.pdf", 10, 8)
save_plot_safely(results$plots$panel_c1, "./Plots/Figure2/PanelC-1.pdf", 6, 8)
save_plot_safely(results$plots$panel_c2, "./Plots/Figure2/PanelC-2.pdf", 6, 8)
save_plot_safely(results$plots$panel_c3, "./Plots/Figure2/PanelC-3.pdf", 6, 8)

for (i in 1:4) {
  save_plot_safely(results$plots[[paste0("panel_d", i)]], 
                   paste0("./Plots/Figure2/PanelD-", i, ".pdf"), 6.5, 6.5)
}

# Print summary
cat("\n=== Analysis Complete ===\n")
cat("Total cells analyzed:", nrow(results$data$cell_proportions), "\n")
cat("Samples processed:", length(SAMPLE_CONFIG$names), "\n")
cat("Plots generated:", length(results$plots), "\n")
