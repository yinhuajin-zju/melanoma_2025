# =============================================================================
# Figure-3 Cell Type Annotation and Visualization
# Author: Hanzhang
# Description: Visualization of data with focus on landscape visualization and
#              cell type distribution of integrated method 3 datasets
# =============================================================================

library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(Nebulosa)

scRNA <- melanoma_seu
celltype_column <- "cell_anno"
DefaultAssay(scRNA) <- "SCT"

# =============================================================================
# Dot Plot Visualization
# =============================================================================

genes_dot <- c("GDF15", "TGFBR2")
genes_dot <- genes_dot[genes_dot %in% rownames(scRNA)]

cell_order <- c("Melanoma cell", "Inflammatory Macrophage", "Macrophage", "Differentiated_KC",
                "Endothelial cell", "Pericyte", "NK", "Treg", "Th", "Tc", "B cell")
scRNA[[celltype_column]] <- factor(scRNA[[celltype_column, drop = TRUE]], levels = rev(cell_order))

p_dot <- DotPlot(scRNA, features = genes_dot, group.by = celltype_column, assay = "SCT") +
  scale_color_gradient2(low = "#0077BB", mid = "white", high = "#CC3311", midpoint = 0, name = "Scaled\nExpression") +
  guides(size = guide_legend(title = 'Percentage\nExpressed')) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  xlab("Genes") + ylab("Cell Types")

print(p_dot)

# =============================================================================
# UMAP Cell Type Visualization
# =============================================================================

plot_umap_celltype <- function(seurat_obj, reduction = "umap", celltype_col = "cell_anno", filter_type = NULL) {
  umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction))
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  umap_coords$celltype <- seurat_obj[[celltype_col, drop = TRUE]]
  umap_coords <- umap_coords[!is.na(umap_coords$celltype), ]
  
  if (!is.null(filter_type)) {
    umap_coords <- if (filter_type == "only_mela") {
      umap_coords[umap_coords$celltype == "Melanoma cell", ]
    } else {
      umap_coords[umap_coords$celltype != "Melanoma cell", ]
    }
  }
  
  if (!nrow(umap_coords)) stop("No cells remain after filtering.")
  
  current_types <- sort(unique(umap_coords$celltype))
  palette <- c("#073A61", "#3F91B4", "#61B673", "#D2DC38", "#C63430", "#f36c24", "#F7B6D3", 
               "#935fa7", "#a54922", "#7F7F7F", "#D4242A", "#17BECF", "#FFBB78", "#98DF8A", 
               "#C5B0D5", "#BCBD22")[seq_along(current_types)]
  
  title_map <- c(
    NULL = "Cell Type Annotation",
    only_mela = "Melanoma Cells Only",
    no_mela = "Non-Melanoma Cells"
  )
  
  ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
    geom_point(alpha = 0.7, shape = 16, size = 0.6) +
    scale_color_manual(values = setNames(palette, current_types)) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      axis.ticks = element_line(color = "black", size = 0.3)
    ) +
    labs(title = title_map[[filter_type]], x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 6)))
}

p_all <- plot_umap_celltype(scRNA, reduction = "umap", celltype_col = "cell_anno")
p_mela <- plot_umap_celltype(scRNA, reduction = "umap", celltype_col = "cell_anno", filter_type = "only_mela")
p_non_m <- plot_umap_celltype(scRNA, reduction = "umap", celltype_col = "cell_anno", filter_type = "no_mela")

ggsave("./Plots/Figure3/Umap_with_mela.png", plot = p_all, width = 10, height = 8, dpi = 300)
ggsave("./Plots/Figure3/Umap_mela_only.png", plot = p_mela, width = 10, height = 8, dpi = 300)
ggsave("./Plots/Figure3/Umap_without_mela.png", plot = p_non_m, width = 10, height = 8, dpi = 300)

# =============================================================================
# Feature Plot Visualization with Density
# =============================================================================

BiocManager::install("Nebulosa")

DefaultAssay(seu) <- "SCT"

p_f1 <- plot_density(seu, features = "CD4", reduction = "umap")
p_f2 <- plot_density(seu, features = "CD22", reduction = "umap")
p_f3 <- plot_density(seu, features = "KLRD1", reduction = "umap")

ggsave(filename = "./Plots/Figure3/F1.pdf", dpi = 100, width = 8, height = 8, plot = p_f1)
ggsave(filename = "./Plots/Figure3/F2.pdf", dpi = 100, width = 8, height = 8, plot = p_f2)
ggsave(filename = "./Plots/Figure3/F3.pdf", dpi = 100, width = 8, height = 8, plot = p_f3)
