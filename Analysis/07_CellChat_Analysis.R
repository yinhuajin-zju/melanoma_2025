# =============================================================================
# Cell-Cell Communication Analysis using CellChat
# Author: Hanzhang
# Description:interaction network construction, pathway analysis, and 
#              visualization of signaling patterns
# =============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(Seurat)
  library(tidyverse)
  library(future)
})

# =============================================================================
# Configuration and Setup
# =============================================================================

plan(sequential)

# =============================================================================
# Data Preparation and CellChat Object Creation
# =============================================================================

expression_matrix <- GetAssayData(seu, assay = "SCT", slot = "data")

# Uncomment if cell annotation needs to be updated
# seu$cell_anno <- seu$majority_voting

cellchat <- createCellChat(
  object = expression_matrix,
  meta = seu@meta.data,
  group.by = "cell_anno"
)

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

cat(paste("Cell groups:", length(unique(cellchat@idents)), "\n"))
cat(paste("Total cells:", ncol(cellchat@data.raw), "\n"))

# =============================================================================
# Data Preprocessing
# =============================================================================

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- smoothData(
  cellchat, 
  adj = PPI.human, 
  alpha = 0.5, 
  normalizeAdjMatrix = "rows"
)

# =============================================================================
# Communication Probability Analysis
# =============================================================================

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# =============================================================================
# Network Analysis and Clustering
# =============================================================================

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# =============================================================================
# Network Visualization
# =============================================================================

netVisual_circle(
  cellchat@net$count, 
  vertex.weight = rowSums(cellchat@net$count),
  weight.scale = TRUE, 
  label.edge = FALSE,
  title.name = "Number of Interactions"
)

netVisual_heatmap(
  cellchat, 
  measure = "weight", 
  font.size = 12, 
  font.size.title = 14
)

netVisual_aggregate(
  object = cellchat, 
  signaling = "GDF", 
  layout = "circle",
  pt.title = 16, 
  vertex.label.cex = 1.2
)

netAnalysis_signalingRole_heatmap(
  cellchat, 
  pattern = "all", 
  width = 20, 
  height = 18, 
  font.size = 10
)

# =============================================================================
# Pathway-Specific Analysis
# =============================================================================

pathways_of_interest <- c("CD99", "APP", "EPHA", "PTPRM", "GDF", "PARs", "CADM")

for (pathway in pathways_of_interest) {
  cat(paste("Visualizing", pathway, "pathway...\n"))
  
  tryCatch({
    netVisual_aggregate(
      cellchat, 
      signaling = pathway, 
      layout = "circle"
    )
  }, error = function(e) {
    cat(paste("Warning: Could not visualize pathway", pathway, ":", e$message, "\n"))
  })
}

# =============================================================================
# Export Results
# =============================================================================

output_dir <- "cellchat_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(
  subsetCommunication(cellchat), 
  file = file.path(output_dir, "CellChat_Interactions.csv"),
  row.names = FALSE
)

write.csv(
  cellchat@netP$pathways, 
  file = file.path(output_dir, "CellChat_Pathways.csv"),
  row.names = FALSE
)

write.csv(
  cellchat@netP$centrality, 
  file = file.path(output_dir, "CellChat_Centrality.csv"),
  row.names = FALSE
)

# =============================================================================
# Gene Expression Visualization
# =============================================================================

GDF_ligand_receptor <- subset(CellChatDB.human$interaction, pathway_name == "GDF")
GDF_genes <- unique(c(GDF_ligand_receptor$ligand, GDF_ligand_receptor$receptor))

cat(paste("GDF pathway genes identified:", length(GDF_genes), "\n"))


gdf_expression_plot <- FeaturePlot(
  seu, 
  features = GDF_genes, 
  reduction = "umap"
) +
  plot_annotation(title = "GDF Pathway Gene Expression")

print(gdf_expression_plot)


gdf_specific_plot <- FeaturePlot(
  seu, 
  features = c("GDF15", "TGFBR2"), 
  reduction = "umap"
) +
  plot_annotation(title = "Key GDF Signaling Components")

print(gdf_specific_plot)

# =============================================================================
# Save CellChat Object
# =============================================================================

saveRDS(cellchat, file = file.path(output_dir, "cellchat_object_annotated.rds"))
