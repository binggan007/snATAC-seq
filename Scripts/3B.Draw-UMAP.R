library(ArchR)
library(gridExtra)
args <- commandArgs(trailingOnly=TRUE)
Accession <- args[2]
proj <- readRDS(args[1])
palette <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"
             ,"#ffff99")
# manual labels store in Labels
# just change the umap color scheme to Labels
p <- getEmbedding(ArchRProj = proj, embedding ="UMAP_Tile")
png(paste0("Unlabeled-UMAP.",Accession,".png"),width=10,height=8,units="in",res=300)
data <- data.frame(table(proj$Clusters_Tile))
names(data) <- c("cellType","count")
barplot <- ggplot(data=data,aes(x=cellType,y=count,fill=cellType))+
           geom_bar(stat="identity")+
           scale_fill_manual(values=palette)+  
           geom_text(aes(label=count),vjust=0.5,hjust=0.5)+
           theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
           theme(legend.position="none")+
           ggtitle(args[2])+
           coord_flip()
umap <- ggplot(p,aes(x=p[,1], y=p[,2], color=proj$Clusters_Tile))+
        geom_point()+
        scale_color_manual(name="Cluster Identity",values=palette)+
        xlab("UMAP 1")+
        ylab("UMAP 2")
grid.arrange(barplot,umap,ncol=1,nrow=2,heights=c(1,3))
dev.off()
print(data)
