library(ArchR)

args <- commandArgs(trailingOnly = TRUE)
proj <- readRDS(args[1])
markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters_Tile", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerGenes <- read.table(file=args[2],sep="\t",header=FALSE)[,1]
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  labelMarkers = markerGenes, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1.25", 
  transpose = TRUE
)

png(paste0("Marker-Heatmap.",args[3],".png"), width=15, height=8, units="in", res=300)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
