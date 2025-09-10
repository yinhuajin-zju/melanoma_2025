# =============================================================================
# ATAC Analysis: Peak-to-Gene Links and Motif Analysis
# Author: Hanzhang
# Description: peak-to-gene linkage, motif enrichment, and gene activity scoring
# =============================================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TFBSTools)
  library(motifmatchr)
  library(JASPAR2022)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# =============================================================================
# Global Configuration
# =============================================================================

set.seed(123)
options(future.globals.maxSize = 20 * 1024^3)
DefaultAssay(seu) <- "ATAC"

# =============================================================================
# Basic Quality Control and Visualization
# =============================================================================

cat("ATAC peaks:", nrow(seu[["ATAC"]]), "\n")
cat("Cells:", ncol(seu[["ATAC"]]), "\n")


qc_plot <- VlnPlot(
  seu, 
  features = c("nCount_ATAC", "nFeature_ATAC"),
  pt.size = 0, 
  ncol = 2, 
  group.by = "cell_anno"
) +
  plot_annotation(title = "ATAC-seq Quality Control Metrics")


umap_plot <- DimPlot(
  seu, 
  reduction = "umap_atac", 
  group.by = "cell_anno", 
  label = TRUE, 
  repel = TRUE
) +
  ggtitle("ATAC-seq UMAP") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(qc_plot)
print(umap_plot)

# =============================================================================
# Peak-to-Gene Linkage Analysis
# =============================================================================

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(seu) <- annotations

seu <- RegionStats(seu, genome = BSgenome.Hsapiens.UCSC.hg38, sep = c("-", "-"))

genes_of_interest <- GDF_genes_raw
seu <- LinkPeaks(
  seu, 
  peak.assay = "ATAC", 
  expression.assay = "SCT", 
  genes.use = genes_of_interest
)

cat("Peak-gene links identified:", length(Links(seu)), "\n")


coverage_gdf15 <- CoveragePlot(
  seu, 
  region = "GDF15", 
  features = "GDF15", 
  expression.assay = "SCT",
  extend.upstream = 5000, 
  extend.downstream = 5000, 
  group.by = "cell_anno"
) +
  ggtitle("GDF15 Coverage Plot") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

coverage_tgfbr2 <- CoveragePlot(
  seu, 
  region = "TGFBR2", 
  features = "TGFBR2", 
  expression.assay = "SCT",
  extend.upstream = 5000, 
  extend.downstream = 5000, 
  group.by = "cell_anno"
) +
  ggtitle("TGFBR2 Coverage Plot") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(coverage_gdf15)
print(coverage_tgfbr2)

# =============================================================================
# Motif Enrichment Analysis
# =============================================================================

pfm <- getMatrixSet(
  JASPAR2022, 
  opts = list(collection = "CORE", tax_group = "vertebrates")
)


seu <- AddMotifs(seu, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
peak_names <- rownames(seu[["ATAC"]])

# TGFBR2 region analysis (chr3:30,640,000-30,670,000)
cat("Analyzing TGFBR2 region motifs...\n")
tgfbr2_region <- peak_names[grepl("^chr3", peak_names)]
tgfbr2_region <- tgfbr2_region[sapply(strsplit(tgfbr2_region, "[-:]"), function(x) {
  coords <- as.numeric(x[2:3])
  length(coords) == 2 && coords[1] >= 30640000 && coords[2] <= 30670000
})]

if (length(tgfbr2_region) > 0) {
  enriched_motifs_tgfbr2 <- FindMotifs(seu, features = head(tgfbr2_region, 10))
  cat("TGFBR2 region enriched motifs:\n")
  print(head(enriched_motifs_tgfbr2))
}

# GDF15 region analysis (chr19:18,370,000-18,390,000)
cat("Analyzing GDF15 region motifs...\n")
gdf15_region <- peak_names[grepl("^chr19", peak_names)]
gdf15_region <- gdf15_region[sapply(strsplit(gdf15_region, "[-:]"), function(x) {
  coords <- as.numeric(x[2:3])
  length(coords) == 2 && coords[1] <= 18390000 && coords[2] >= 18370000
})]

if (length(gdf15_region) > 0) {
  enriched_motifs_gdf15 <- FindMotifs(seu, features = gdf15_region, assay = "ATAC")
  cat("GDF15 region enriched motifs:\n")
  print(head(enriched_motifs_gdf15))
}


motif_plot_tgfbr2 <- MotifPlot(
  seu, 
  motifs = head(rownames(enriched_motifs_tgfbr2), 6), 
  ncol = 3
) +
  plot_annotation(title = "TGFBR2 Region Enriched Motifs")

motif_plot_gdf15 <- MotifPlot(
  seu, 
  motifs = head(rownames(enriched_motifs_gdf15), 8), 
  ncol = 4
) +
  plot_annotation(title = "GDF15 Region Enriched Motifs")

print(motif_plot_tgfbr2)
print(motif_plot_gdf15)

# =============================================================================
# Transcription Factor Footprinting Analysis
# =============================================================================


footprint_motifs <- c("TFAP4", "MSC", "MYF6")


seu <- Footprint(
  seu, 
  motif.name = footprint_motifs, 
  genome = BSgenome.Hsapiens.UCSC.hg38
)


footprint_plot <- PlotFootprint(seu, features = footprint_motifs) +
  plot_annotation(title = "Transcription Factor Footprints")

print(footprint_plot)

# =============================================================================
# Gene Activity Score Analysis
# =============================================================================

gene_activities <- GeneActivity(seu)


seu[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)
seu <- NormalizeData(
  seu, 
  assay = "ACTIVITY", 
  normalization.method = "LogNormalize",
  scale.factor = median(seu$nCount_ACTIVITY)
)

genes_for_activity <- c("TFAP4", "MSC", "MYF6")

DefaultAssay(seu) <- "ACTIVITY"
activity_plot <- FeaturePlot(
  seu, 
  features = genes_for_activity, 
  pt.size = 0.1, 
  max.cutoff = "q95", 
  ncol = 3, 
  reduction = "umap"
) +
  plot_annotation(title = "Gene Activity Scores")

DefaultAssay(seu) <- "SCT"
expression_plot <- FeaturePlot(
  seu, 
  features = genes_for_activity, 
  pt.size = 0.1, 
  max.cutoff = "q95", 
  ncol = 3, 
  reduction = "umap"
) +
  plot_annotation(title = "Gene Expression")

print(activity_plot)
print(expression_plot)

# =============================================================================
# Save Results
# =============================================================================

saveRDS(seu, "M3_with_atac_analysis_final.rds")
