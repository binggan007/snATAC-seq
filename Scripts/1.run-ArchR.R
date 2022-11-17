# Yu Fu
# Weng Lab
# This scripts takes fragment file as input and creates an ArchR project with the following parameters(other parameters as default):
        # bsgenome with v29 annotation and ENCODE blacklist regions
        # TSS threshold is 7
        # UMAP parameters: minDist=0.1 and nNeighbors=10
        # Embedding with two resolution: UMAP_Tile(resolution=2) and UMAP_Tile_Low(resolution=.2)

library(ArchR)
library(GenomicRanges)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
out <- args[2]
names <- args[3]

# our workflow with tile matrix; w/ bs gene annotation
blacklist_path <- "ENCFF356LFX.bed"
regions <- read.table(blacklist_path, sep = '\t', header = FALSE)
colnames(regions) <- c('chr','start','end')
blacklist <- GRanges(regions)

bsgenome_name="BSgenome.Hsapiens.UCSC.hg38"
library(bsgenome_name, character.only = TRUE)
bsgenome <- get(bsgenome_name)

chromSizes <- GRanges(names(seqlengths(bsgenome)), IRanges(1, seqlengths(bsgenome)))
chromSizes <- filterChrGR(chromSizes, remove = c("chrM"))
seqlengths(chromSizes) <- end(chromSizes)

genome_annotation <- createGenomeAnnotation(
        genome = bsgenome,
        chromSizes = chromSizes,
        blacklist = blacklist
    )

library(GenomicFeatures)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
gene_annotation <- createGeneAnnotation(TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,OrgDb=org.Hs.eg.db)
ArrowFiles = createArrowFiles(inputFiles = input, sampleNames = names,
			       genomeAnnotation = genome_annotation,
			       geneAnnotation = gene_annotation,
			       force=TRUE,addTileMat = TRUE, addGeneScoreMat = TRUE,minTSS = 7, minFrags = 1000)

proj = ArchRProject(ArrowFiles = ArrowFiles,genomeAnnotation = genome_annotation, geneAnnotation = gene_annotation,  outputDirectory = out, copyArrows = FALSE)
proj = addDoubletScores(proj,k=10,knnMethod="UMAP",LSIMethod=1)
proj = filterDoublets(proj)

# IterativeLSI resolution = 2
proj = addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI_Tile",iterations = 2, clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30, force=TRUE)
proj = addClusters(input = proj, reducedDims = "IterativeLSI_Tile",method = "Seurat", name = "Clusters_Tile", resolution = 0.8, force=TRUE)
proj = addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI_Tile", name = "UMAP_Tile", nNeighbors = 10, minDist = 0.1, metric = "cosine", force=TRUE)
getGroupBW(ArchRProj=proj, groupBy="Clusters_Tile")

# IterativeLSI resolution = .2
#proj = addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI_Tile_Low",iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30, force=TRUE)
#proj = addClusters(input = proj, reducedDims = "IterativeLSI_Tile_Low",method = "Seurat", name = "Clusters_Tile_Low", resolution = 0.8, force=TRUE)
#proj = addUMAP(ArchRProj  = proj, reducedDims = "IterativeLSI_Tile_Low",name = "UMAP_Tile_Low",nNeighbors = 10, minDist = 0.1, metric = "cosine", force=TRUE)
#getGroupBW(ArchRProj=proj, groupBy="Clusters_Tile_Low")

saveArchRProject(ArchRProj = proj, outputDirectory = out, overwrite = TRUE, load = FALSE)
