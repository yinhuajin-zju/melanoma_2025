# =============================================================================
# Figure-4 Featureplots of interested genes
# Author: Hanzhang
# Description: Modification of featureplots
# =============================================================================

#############################################
# Figure 4: FeaturePlots modification
#############################################
library(Seurat)
library(ggplot2)
library(scCustomize)

# 状态：输入 scRNA 对象
my_colors <- colorRampPalette(c("#ED213A", "#93291E"))(50)
point_size <- 0.5

genes_feat <- c("GDF15", "TGFBR2", "MYF6", "TFAP4", "MSC", "TCF3")
plots_feat <- list()

for (g in genes_feat) {
  if (g %in% rownames(scRNA)) {
    p <- FeaturePlot_scCustom(
      seurat_object = scRNA,
      colors_use    = my_colors,
      features      = g,
      pt.size       = point_size
    ) + labs(x = "UMAP1", y = "UMAP2") +
      theme(axis.title = element_text(face = "bold", size = 12))
    plots_feat[[g]] <- p
    print(p)
  }
}

dir.create("./Plots/Figure4", recursive = TRUE, showWarnings = FALSE)
for (g in names(plots_feat)) {
  ggsave(sprintf("./Plots/Figure4/FeaturePlot_%s.png", g),
         plot = plots_feat[[g]], width = 10, height = 8, dpi = 300)
}
