library(ArchR)
library(viridis)
library(tidyverse)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly=TRUE)
proj <- loadArchRProject(path=args[1])
markers <- scan(args[2],what="")
markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix",groupBy = "Clusters_Tile", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS,cutOff = "FDR <= 0.1 & Log2FC >= 1.25")
seGS <- getGroupSE(ArchRProj=proj, useMatrix = "GeneScoreMatrix", groupBy = "Clusters_Tile")
groupGeneScore <- assay(seGS)
rownames(groupGeneScore) <- rowData(seGS)$name
markersGeneScore <- groupGeneScore[rownames(groupGeneScore) %in% markers,]
mode(markersGeneScore) <- 'numeric'
markersGeneScore <- as.data.frame(markersGeneScore)
col_fun = circlize::colorRamp2(c(0,0.3,1.5), viridis(20)[c(1,10, 20)])
cell_fun = function(j, i, x, y, w, h, fill){
	grid.rect(x = x, y = y, width = w, height = h,gp = gpar(col = NA, fill = NA))
	grid.circle(x=x,y=y,r= unit(2, "mm"),gp = gpar(fill = col_fun(markersGeneScore[i, j]), col = NA))}
png(paste0("Dotplot.",args[3],".png"),width=10, height=8, units="in", res=300)
Heatmap(markersGeneScore,
	heatmap_legend_param=list(title="expression"),
	column_title = "Clustered Dotplot", 
	col=col_fun,
	rect_gp = gpar(type = "none"),
	cell_fun = cell_fun,
	row_names_gp = gpar(fontsize = 10),
	row_km = 5,
	border = "black")
