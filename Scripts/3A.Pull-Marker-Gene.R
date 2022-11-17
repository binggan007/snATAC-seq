#This scripts will pull out the marker features in each cluster by FDR<=0.2 and Log2FC>=1.25 threshold and save the table.
#Cluster with resolution=2
library(parallel)
library(ArchR)

args <- commandArgs(trailingOnly = TRUE)
proj <- loadArchRProject(args[1])
markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", 
			         groupBy = "Clusters_Tile", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.1 & Log2FC >= 1.25")
clusters  <- names(markerList)
for (name in clusters){
        C <- markerList[[name]]
	write.table(C,paste0(args[2],"/",name,"-diff-genes.txt"),sep="\t",row.names=FALSE,quote=FALSE)
}
