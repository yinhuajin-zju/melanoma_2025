# =============================================================================
# GDF15-Associated Differential Peak and Motif Analysis
# Author: Hanzhang
# Description: Analysis of differential chromatin accessibility
#              and motif enrichment between GDF15-positive and GDF15-negative
#              melanoma cells
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2022)
})

# =============================================================================
# Data Preparation and Cell Subset Selection
# =============================================================================

# Note: Input 'seu' object must contain both ATAC and SCT assays
if (!all(c("ATAC", "SCT") %in% names(seu@assays))) {
  stop("Input Seurat object must contain both ATAC and SCT assays")
}

# Subset to melanoma cells only
melanoma_seu <- subset(seu, subset = cell_anno == "Melanoma cell")

cat(paste("Melanoma cells selected:", ncol(melanoma_seu), "\n"))

# =============================================================================
# GDF15 Expression-Based Cell Classification
# =============================================================================

DefaultAssay(melanoma_seu) <- "SCT"

gdf15_expression <- GetAssayData(melanoma_seu, assay = "SCT", slot = "data")["GDF15", ]
melanoma_seu$gdf15_status <- ifelse(gdf15_expression > 0, "GDF15_Pos", "GDF15_Neg")

classification_table <- table(melanoma_seu$gdf15_status)
print(classification_table)

cat(paste("GDF15_Pos cells:", classification_table["GDF15_Pos"], "\n"))
cat(paste("GDF15_Neg cells:", classification_table["GDF15_Neg"], "\n"))

# =============================================================================
# Differential Accessibility Analysis
# =============================================================================

DefaultAssay(melanoma_seu) <- "ATAC"
Idents(melanoma_seu) <- "gdf15_status"

differential_peaks <- FindMarkers(
  object = melanoma_seu,
  ident.1 = "GDF15_Pos",
  ident.2 = "GDF15_Neg",
  test.use = "LR",                    # Likelihood ratio test
  latent.vars = "nCount_ATAC",        # Control for sequencing depth
  logfc.threshold = 0.25,             # Minimum log fold change
  min.pct = 0.05                      # Minimum percentage of cells expressing
)

significant_peaks_pos <- differential_peaks %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0)

cat(paste("Significant peaks enriched in GDF15_Pos group:", nrow(significant_peaks_pos), "\n"))
write.csv(
  differential_peaks, 
  "GDF15_differential_peaks_results.csv", 
  row.names = TRUE
)

# =============================================================================
# Motif Enrichment Analysis for GDF15_Pos Peaks
# =============================================================================

if (nrow(significant_peaks_pos) > 0) {
  # Load motif database
  pfm <- getMatrixSet(
    JASPAR2022, 
    opts = list(collection = "CORE", tax_group = "vertebrates")
  )
  
  if (is.null(melanoma_seu[["ATAC"]]@motifs)) {
    cat("Adding motif information to Seurat object...\n")
    melanoma_seu <- AddMotifs(
      melanoma_seu, 
      genome = BSgenome.Hsapiens.UCSC.hg38, 
      pfm = pfm
    )
  }
  
  enriched_motifs <- FindMotifs(
    melanoma_seu, 
    features = rownames(significant_peaks_pos)
  )
  
  write.csv(
    enriched_motifs, 
    "GDF15_Pos_vs_Neg_enriched_motifs.csv", 
    row.names = FALSE
  )
  
  cat("Top enriched motifs:\n")
  print(head(enriched_motifs))
  
} else {
  cat("No significant peaks found for motif analysis\n")
}

# =============================================================================
# Coverage Plot Comparison
# =============================================================================

# Define peak region of interest (GDF15 locus)
peak_region <- "chr19-18385688-18386608"


subset_for_plot <- subset(melanoma_seu, subset = gdf15_status %in% c("GDF15_Pos", "GDF15_Neg"))
subset_for_plot$gdf15_status <- factor(
  subset_for_plot$gdf15_status, 
  levels = c("GDF15_Pos", "GDF15_Neg")
)


coverage_plot <- CoveragePlot(
  object = subset_for_plot,
  region = peak_region,
  extend.upstream = 10000,
  extend.downstream = 10000,
  annotation = TRUE,
  peaks = TRUE,
  group.by = "gdf15_status"
) & 
  scale_fill_manual(
    values = c(
      "GDF15_Pos" = rgb(82, 170, 220, maxColorValue = 255),
      "GDF15_Neg" = rgb(236, 110, 102, maxColorValue = 255)
    )
  ) &
  theme(
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
    legend.position = "none"
  )

print(coverage_plot)
ggsave("GDF15_coverage_comparison.pdf", plot = coverage_plot, width = 12, height = 8)

# =============================================================================
# Single Peak Motif Analysis
# =============================================================================

single_peak_motifs <- FindMotifs(melanoma_seu, features = peak_region)
write.csv(
  single_peak_motifs, 
  paste0("motif_enrichment_", gsub("[-:]", "_", peak_region), ".csv"), 
  row.names = FALSE
)


print(head(single_peak_motifs))

# =============================================================================
# Regional Motif Analysis and Visualization
# =============================================================================


output_directory <- "./Plots/Figure4"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# Define target genomic region (GDF15 locus)
target_region_start <- 18388000
target_region_end <- 18392500


chr19_peaks <- grep("chr19", rownames(melanoma_seu), value = TRUE)
regional_peaks <- chr19_peaks[sapply(chr19_peaks, function(peak_name) {
  coordinates <- as.numeric(strsplit(peak_name, "-")[[1]][2:3])
  coordinates[1] <= target_region_end & coordinates[2] >= target_region_start
})]

cat(paste("Found", length(regional_peaks), "peaks in target region\n"))
regional_motif_results <- FindMotifs(melanoma_seu, features = regional_peaks)


write.csv(
  regional_motif_results, 
  file.path(output_directory, "chr19_regional_motifs.csv"),
  row.names = FALSE
)

# =============================================================================
# Top Motif Visualization and Footprinting
# =============================================================================

top_motifs <- head(regional_motif_results, 6)


for (i in seq_len(nrow(top_motifs))) {
  motif_identifier <- top_motifs$motif[i]
  motif_name_clean <- toupper(top_motifs$motif.name[i])
  

    motif_logo <- MotifPlot(melanoma_seu, motifs = motif_identifier) +
    ggtitle(motif_name_clean) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text = element_blank()
    )
  
  safe_filename <- gsub("[^A-Za-z0-9_\\-]", "_", motif_name_clean)
  ggsave(
    file.path(output_directory, paste0("MotifLogo_", safe_filename, ".pdf")),
    plot = motif_logo, 
    width = 4, 
    height = 3
  )
}

target_motif_names <- top_motifs$motif.name


melanoma_seu <- Footprint(
  melanoma_seu, 
  motif.name = target_motif_names, 
  genome = BSgenome.Hsapiens.UCSC.hg38
)


footprint_plots <- PlotFootprint(melanoma_seu, features = target_motif_names)


for (i in seq_along(footprint_plots)) {
  footprint_name <- toupper(target_motif_names[i])
  
  footprint_plot <- footprint_plots[[i]] + 
    theme_void() +
    ggtitle(footprint_name) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
  

    safe_fp_name <- gsub("::", "_", footprint_name)
  
  ggsave(
    file.path(output_directory, paste0("footprint_", safe_fp_name, ".pdf")),
    footprint_plot, 
    width = 8, 
    height = 6
  )
}

