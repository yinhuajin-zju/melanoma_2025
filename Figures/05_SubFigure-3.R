# =============================================================================
# SubFigure-3: Marker Gene Visualization
# Author: Hanzhang
# Description: Marker Gene Visualization by using Nebulosa
# =============================================================================

library(Seurat)
library(Nebulosa)
library(ggplot2)

# =============================================================================
# Data Loading
# =============================================================================

M1_nofix <- readRDS("./Data/M1_nofix_wnn.RDS")

M3_0045 <- readRDS("./Data/updated_0045_wnn.RDS")
M3_0051 <- readRDS("./Data/updated_0051_wnn.RDS")
M3_0060 <- readRDS("./Data/updated_0060_wnn.RDS")
M3_0067 <- readRDS("./Data/updated_0067_wnn.RDS")

# =============================================================================
# UMAP Cluster Visualization
# =============================================================================

DimPlot(M3_0045, reduction = "wnn.umap")
DimPlot(M3_0051, reduction = "wnn.umap")
DimPlot(M3_0060, reduction = "wnn.umap")
DimPlot(M3_0067, reduction = "wnn.umap")

# =============================================================================
# Marker Gene Density Plots for M3_0045
# =============================================================================

BiocManager::install("Nebulosa")

DefaultAssay(M3_0045) <- "SCT"

p_sf1 <- plot_density(M3_0045, features = "CD4", reduction = "wnn.umap")
p_sf2 <- plot_density(M3_0045, features = "CD22", reduction = "wnn.umap")
p_sf3 <- plot_density(M3_0045, features = "CD68", reduction = "wnn.umap")
p_sf4 <- plot_density(M3_0045, features = "CD14", reduction = "wnn.umap")
p_sf5 <- plot_density(M3_0045, features = "KLRD1", reduction = "wnn.umap")
p_sf6 <- plot_density(M3_0045, features = "NCAM1", reduction = "wnn.umap")
p_sf7 <- plot_density(M3_0045, features = "NCR1", reduction = "wnn.umap")
p_sf8 <- plot_density(M3_0045, features = "KLRK1", reduction = "wnn.umap")

ggsave(filename = "./Plots/SubFigure3/CD4.png", plot = p_sf1, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/CD22.png", plot = p_sf2, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/CD68.png", plot = p_sf3, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/CD14.png", plot = p_sf4, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/KLRD1.png", plot = p_sf5, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/NCAM1.png", plot = p_sf6, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/NCR1.png", plot = p_sf7, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/KLRK1.png", plot = p_sf8, width = 10, height = 8, dpi = 300)

# =============================================================================
# Marker Gene Density Plots for Integrated Data
# =============================================================================

seu <- readRDS("./M3_with_atac_analysis_final.rds")

DefaultAssay(seu) <- "SCT"

p_sf9 <- plot_density(seu, features = "CD68", reduction = "umap")
p_sf10 <- plot_density(seu, features = "CD14", reduction = "umap")
p_sf11 <- plot_density(seu, features = "NCAM1", reduction = "umap")
p_sf12 <- plot_density(seu, features = "NCR1", reduction = "umap")
p_sf13 <- plot_density(seu, features = "KLRK1", reduction = "umap")

ggsave(filename = "./Plots/SubFigure3/CD68_M3.png", plot = p_sf9, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/CD14_M3.png", plot = p_sf10, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/NCAM1_M3.png", plot = p_sf11, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/NCR1_M3.png", plot = p_sf12, width = 10, height = 8, dpi = 300)
ggsave(filename = "./Plots/SubFigure3/KLRK1_M3.png", plot = p_sf13, width = 10, height = 8, dpi = 300)
