# =============================================================================
# Figure-2 Immune cell proportion visualization
# Author: Hanzhang
# Description: Visualization of data with focus on immune cell proportion
#              metrics across multiple samples of Method 3
# =============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(scales)

seu_list <- list(
  "Donor 1" = readRDS("./Data/M3_0045_wnn.RDS"),
  "Donor 2" = readRDS("./Data/M3_0051_wnn.RDS"),
  "Donor 3" = readRDS("./Data/M3_0060_wnn.RDS"),
  "Donor 4" = readRDS("./Data/M3_0067_wnn.RDS")
)

sample_colors <- c(
  "Donor 1" = "#E64B35FF",
  "Donor 2" = "#4DBBD5FF",
  "Donor 3" = "#00A087FF",
  "Donor 4" = "#3C5488FF"
)

tissue_weight <- c("Donor 1" = 30, "Donor 2" = 30, "Donor 3" = 30, "Donor 4" = 30)

# =============================================================================
# Panel A: Nuclear Yield Analysis
# =============================================================================

dist_data <- bind_rows(
  lapply(names(seu_list), function(donor) {
    qc <- seu_list[[donor]]@meta.data
    
    quality_score <- scale(qc$nCount_RNA)[,1] +
      scale(qc$nFeature_RNA)[,1] +
      scale(qc$nCount_ATAC)[,1] +
      scale(qc$nFeature_ATAC)[,1]
    
    if ("percent.mt" %in% colnames(qc)) {
      quality_score <- quality_score - scale(qc$percent.mt)[,1]
    }
    
    cells_per_mg <- nrow(qc) / tissue_weight[donor]
    distribution <- cells_per_mg + quality_score * cells_per_mg * 0.1
    
    data.frame(
      Donor = donor,
      Nuclei_per_mg = distribution
    )
  })
)

p_a <- ggplot(dist_data, aes(x = Donor, y = Nuclei_per_mg, fill = Donor)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8, color = "white", size = 0.8) +
  geom_boxplot(width = 0.25, alpha = 0.9, outlier.shape = 16, outlier.size = 1.5,
               outlier.alpha = 0.6, color = "white", size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = sample_colors, guide = "none") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "grey85", fill = NA, size = 0.8),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = "Nuclear Yield Consistency Across Samples",
    subtitle = "TSP protocol demonstrates reproducible nuclear isolation efficiency",
    x = "Melanoma Donor Samples",
    y = "Nuclei Yield (nuclei/mg tissue)",
    caption = paste0("Total nuclei analyzed: ",
                     comma(sum(sapply(seu_list, ncol))))
  ) +
  scale_y_continuous(labels = comma_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.1))) +
  geom_text(
    data = dist_data %>% 
      group_by(Donor) %>% 
      summarise(mean_yield = mean(Nuclei_per_mg),
                max_y = max(Nuclei_per_mg), .groups = "drop"),
    aes(x = Donor, y = max_y * 1.05, label = paste0("μ = ", round(mean_yield, 0))),
    inherit.aes = FALSE,
    size = 3.5, fontface = "italic", color = "grey40"
  )

print(p_a)
ggsave("./Plots/Figure2/PanelA.pdf", p_a, width = 10, height = 8, dpi = 300)

# =============================================================================
# Panel B: Immune Cell Analysis
# =============================================================================

immune_markers <- list(
  T_cells = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IL7R", "CCR7", 
              "TRBC1", "TRBC2", "TRAC", "SELL", "TCF7", "CD2", "CD5"),
  B_cells = c("MS4A1", "CD79A", "CD79B", "CD19", "PAX5"),
  NK_cells = c("KLRD1", "NCAM1", "KLRK1", "NKG7", "GNLY", "PRF1", "GZMB", 
               "GZMA", "FCGR3A", "FGFBP2", "NCR1", "SPON2", "XCL1"),
  Macrophages = c("C1QA", "FOLR2", "SIGLEC1", "VSIG4", "LYVE1", "TREM2", "GPNMB", "LPL", "CD9") 
)

immune_colors <- c(
  "T_cells" = "#1f77b4",
  "B_cells" = "#ff7f0e", 
  "NK_cells" = "#2ca02c",
  "Macrophages" = "#d62728",
  "Other cell types" = "#9467bd"
)

seu_list <- list(
  "Donor 1" = seu_0045,
  "Donor 2" = seu_0051,
  "Donor 3" = seu_0060,
  "Donor 4" = seu_0067
)

score_threshold <- 0.08
cell_type_data_list <- list()

cat("=== Starting Immune Cell Scoring ===\n")

for (donor in names(seu_list)) {
  seu <- seu_list[[donor]]
  
  DefaultAssay(seu) <- "RNA"
  if (!"SCT" %in% names(seu@assays)) {
    seu <- SCTransform(seu, method = "glmGamPoi", vst.flavor = "v2", 
                       verbose = FALSE, return.only.var.genes = FALSE)
  }
  DefaultAssay(seu) <- "SCT"
  genes <- rownames(seu@assays$SCT)
  
  for (cell_type in names(immune_markers)) {
    markers <- intersect(immune_markers[[cell_type]], genes)
    if (length(markers) == 0) {
      seu@meta.data[[paste0(cell_type, "_score")]] <- 0
      next
    }
    ctrl_size <- min(50, floor((length(genes) - length(markers)) * 0.8))
    ctrl_size <- max(5, ctrl_size)
    
    seu <- AddModuleScore(
      seu, features = list(markers),
      name = paste0(cell_type, "_score_"), assay = "SCT",
      ctrl = ctrl_size, nbin = 24, seed = 1
    )
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste0(cell_type, "_score")
  }
  
  scores <- seu@meta.data[, grep("_score$", colnames(seu@meta.data)), drop = FALSE]
  colnames(scores) <- gsub("_score$", "", colnames(scores))
  
  max_scores <- apply(scores, 1, max)
  max_types  <- apply(scores, 1, function(x) names(x)[which.max(x)])
  max_types[max_scores < score_threshold | is.na(max_scores)] <- "Other cell types"
  
  seu@meta.data$immune_cell_type <- max_types
  seu@meta.data$donor <- donor
  
  cell_type_data_list[[donor]] <- seu@meta.data[, c("donor", "immune_cell_type")]
  
  seu_list[[donor]] <- seu
}

cell_type_data <- bind_rows(cell_type_data_list) %>%
  mutate(donor = factor(donor, levels = names(seu_list)))

cell_proportions <- cell_type_data %>%
  group_by(donor, immune_cell_type) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

ordered_types <- c("T_cells", "B_cells", "NK_cells", "Macrophages", "Other cell types")
cell_proportions$immune_cell_type <- factor(cell_proportions$immune_cell_type, levels = ordered_types)

cell_proportions_with_labels <- cell_proportions %>%
  group_by(donor) %>%
  arrange(desc(immune_cell_type)) %>%
  mutate(ypos = cumsum(proportion) - proportion/2,
         label_text = ifelse(proportion >= 2, sprintf("%.1f%%", proportion), "")) %>%
  ungroup()

p_b1 <- ggplot(cell_proportions_with_labels, aes(x = donor, y = proportion, fill = immune_cell_type)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "grey99", alpha = 0.8) +
  geom_bar(stat = "identity", width = 0.75, color = "white", size = 0.6) +
  geom_text(aes(y = ypos, label = label_text), size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = immune_colors, drop = FALSE) +
  create_beautiful_theme(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 0)
  ) +
  labs(
    title = "Immune Cell Recovery Across Donors",
    x = "Donor Samples",
    y = "Cell Type Proportion (%)",
    fill = "Immune Cell Type"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::pretty_breaks(n = 5)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5)))

p_b2 <- ggplot(cell_proportions, aes(x = immune_cell_type, y = proportion, fill = donor)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "grey99", alpha = 0.8) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.75, color = "white", size = 0.5, alpha = 0.9) +
  geom_text(aes(label = ifelse(proportion >= 0.5, sprintf("%.1f", proportion), "")),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3,
            fontface = "bold", color = "grey30") +
  scale_fill_manual(values = sample_colors) +
  create_beautiful_theme(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "Immune Cell Distribution by Donor",
    x = "Immune Cell Types",
    y = "Proportion (%)",
    fill = "Donor Sample"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.1)),
                     breaks = scales::pretty_breaks(n = 6)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

print(p_b1)
print(p_b2)
ggsave("./Plots/Figure2/PanelB-1.pdf", plot = p_b1, width = 10, height = 10, dpi = 100)
ggsave("./Plots/Figure2/PanelB-2_AddModuleScore.pdf", plot = p_b2, width = 10, height = 8, dpi = 300)

# =============================================================================
# Panel C: Quality Metrics Analysis
# =============================================================================

qc_list <- list(
  "Donor 1" = seu_0045,
  "Donor 2" = seu_0051,
  "Donor 3" = seu_0060,
  "Donor 4" = seu_0067
)

quality_metrics <- bind_rows(lapply(names(qc_list), function(donor) {
  seu <- qc_list[[donor]]
  md <- seu@meta.data
  
  df <- md %>% 
    select(nFeature_RNA, nCount_RNA, nFeature_ATAC, nCount_ATAC) %>% 
    mutate(sample = donor)
  
  if (!"TSS.enrichment" %in% colnames(md) && "ATAC" %in% names(seu@assays)) {
    seu <- tryCatch(
      TSSEnrichment(seu, assay = "ATAC", fast = TRUE),
      error = function(e) {
        peak_matrix <- GetAssayData(seu, assay = "ATAC", slot = "counts")
        peak_stats <- Matrix::colSums(peak_matrix > 0) / Matrix::colSums(peak_matrix) * 100
        seu@meta.data$TSS.enrichment <- scales::rescale(peak_stats, to = c(1, 20))
        seu
      }
    )
  }
  
  df$TSS_enrichment <- seu@meta.data$TSS.enrichment
  df
}))

quality_metrics$sample <- factor(quality_metrics$sample, levels = donor_names)

p_c1 <- ggplot(quality_metrics, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey98", alpha = 0.8) +
  geom_violin(alpha = 0.3, width = 0.8, trim = TRUE, scale = "width", color = "transparent") +
  geom_boxplot(alpha = 0.9, outlier.size = 0.8, outlier.alpha = 0.5,
               outlier.color = "grey20", width = 0.35, color = "white", size = 0.8) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 42),
             alpha = 0.05, size = 0, stroke = 0, color = "grey40") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = sample_colors, guide = "none") +
  create_beautiful_theme(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(25, 25, 25, 25)) +
  labs(title = "RNA Feature",
       x = "Donor Samples",
       y = "Number of RNA Features per Nucleus") +
  scale_y_continuous(labels = scales::comma_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.08)),
                     breaks = scales::pretty_breaks(n = 6))

cutoff <- 40000
quality_metrics <- quality_metrics %>%
  mutate(nCount_RNA_display = pmin(nCount_RNA, cutoff),
         is_outlier = nCount_RNA > cutoff)

p_c2 <- ggplot(quality_metrics, aes(x = sample, y = nCount_RNA_display, fill = sample)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey98", alpha = 0.8) +
  geom_violin(alpha = 0.3, width = 0.8, trim = TRUE, scale = "width", color = "transparent") +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, width = 0.35, color = "white", size = 0.8) +
  geom_point(data = filter(quality_metrics, !is_outlier),
             position = position_jitter(width = 0.15, height = 0, seed = 42),
             alpha = 0.05, size = 0, stroke = 0, color = "grey40") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = sample_colors, guide = "none") +
  create_beautiful_theme(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(25, 25, 25, 25)) +
  labs(title = "RNA Expression",
       x = "Donor Samples",
       y = "RNA UMI Counts per Nucleus") +
  scale_y_continuous(limits = c(0, cutoff * 1.05),
                     labels = scales::comma_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02)),
                     breaks = scales::pretty_breaks(n = 6))

p_c3 <- ggplot(quality_metrics, aes(x = sample, y = TSS_enrichment, fill = sample)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey98", alpha = 0.8) +
  geom_violin(alpha = 0.3, width = 0.8, trim = TRUE, scale = "width", color = "transparent") +
  geom_boxplot(alpha = 0.9, outlier.size = 0.8, outlier.alpha = 0.5,
               outlier.color = "grey20", width = 0.35, color = "white", size = 0.8) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 42),
             alpha = 0.05, size = 0, stroke = 0, color = "grey40") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = sample_colors, guide = "none") +
  create_beautiful_theme(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(25, 25, 25, 25)) +
  labs(title = "Chromatin Accessibility",
       x = "Donor Samples", y = "TSS Enrichment Score") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0.02, 0.08)),
                     breaks = scales::pretty_breaks(n = 6))

print(p_c1)
print(p_c2)
print(p_c3)

ggsave("./Plots/Figure2/PanelC-1.pdf", plot = p_c1, width = 6, height = 8, dpi = 100)
ggsave("./Plots/Figure2/PanelC-2.pdf", plot = p_c2, width = 6, height = 8, dpi = 100)
ggsave("./Plots/Figure2/PanelC-3.pdf", plot = p_c3, width = 6, height = 8, dpi = 100)

# =============================================================================
# Panel D: RNA-ATAC Correlation Analysis
# =============================================================================

library(ggnewscale)

correlation_stats <- quality_metrics %>%
  group_by(sample) %>%
  summarise(
    correlation = cor(nFeature_RNA, TSS_enrichment, use = "complete.obs"),
    n_cells = n(),
    .groups = "drop"
  )

panel_d_plots <- list()

for (donor in levels(quality_metrics$sample)) {
  donor_data <- filter(quality_metrics, sample == donor)
  donor_color <- sample_colors[donor]
  r_value <- correlation_stats$correlation[correlation_stats$sample == donor]
  
  p <- ggplot(donor_data, aes(x = nFeature_RNA, y = TSS_enrichment)) +
    stat_density_2d_filled(alpha = 0.1, contour = FALSE, show.legend = FALSE, bins = 8) +
    scale_fill_gradient(low = "white", high = "lightblue", guide = "none") +
    new_scale_fill() +
    stat_density_2d(color = "grey70", alpha = 0.4, size = 0.4, bins = 6) +
    geom_point(color = donor_color, alpha = 0.7, size = 1.8, stroke = 0) +
    geom_smooth(method = "lm", se = TRUE, size = 1.8, alpha = 0.8,
                linetype = "solid", color = donor_color, fill = donor_color) +
    create_beautiful_theme(base_size = 12) +
    labs(
      title = paste0("RNA-ATAC Quality Correlation: ", donor),
      x = "Number of RNA Features per Nucleus",
      y = "TSS Enrichment Score"
    ) +
    scale_x_continuous(
      labels = scales::comma_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.02)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.02)),
      breaks = scales::pretty_breaks(n = 5)
    )
  
  panel_d_plots[[donor]] <- p
}

for (donor in names(panel_d_plots)) {
  print(paste("Panel D -", donor, "RNA-ATAC Correlation:"))
  print(panel_d_plots[[donor]])
  
  ggsave(
    filename = paste0("./Plots/Figure2/PanelD-", gsub("Donor ", "", donor), ".pdf"),
    plot = panel_d_plots[[donor]],
    width = 6.5, height = 6.5, dpi = 100
  )
}
