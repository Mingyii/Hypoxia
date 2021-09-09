setwd('/data/5.T cell trajectory')
library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle)
source("../4.Cell interaction/scores_RNA.r")
HypoxiaMarkers <- read.delim("../1.tumorigenesis/HypoxiaMarkerSignature.txt",header=F)$V1
human_mouseTrans <- read.csv('mouse.human.homology.genes.csv')
HypoxiaMarkersM <- human_mouseTrans[match(HypoxiaMarkers,human_mouseTrans$hgnc_symbol),'mgi_symbol']
HypoxiaMarkersML <- list(HypoxiaMarkers=HypoxiaMarkersM)

##load mouse T cell Seurat
##data from: Interpretation of T cell states from single-cell transcriptomics data using reference atlases
Tcell <- readRDS('TILAtlas_mouse.rds')
Idents(Tcell) <- 'functional.cluster'
CD8Tcell <- subset(Tcell,idents  = c('CD8_NaiveLike','CD8_Tex','CD8_Tpex','CD8_EffectorMemory','CD8_EarlyActiv')) #subset CD8T cells

pdf('3.figure/6.NCmouse_CD8Tcell_dimplot.pdf',5,3.5)
DimPlot(CD8Tcell)
dev.off()

HPscores <- score_cells(seur=CD8Tcell, names=HypoxiaMarkersML, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
CD8Tcell$hypoxia <- HPscores[match(rownames(CD8Tcell@meta.data),rownames(HPscores)),'HypoxiaMarkers']
CD8Tcell_meta <- CD8Tcell@meta.data

CD8Tcell_meta$functional.cluster <- as.character(CD8Tcell_meta$functional.cluster)
x_order <- lapply(split(CD8Tcell_meta$hypoxia,CD8Tcell_meta$functional.cluster), median) %>% unlist()
x_order <- names(x_order)[order(x_order)]


p <- ggplot(CD8Tcell_meta,aes(x=functional.cluster,y=hypoxia))+
  geom_boxplot(aes(fill=functional.cluster),width = 0.6)+
  scale_x_discrete(limits=x_order,label=paste0(x_order,"(",table(CD8Tcell_meta$functional.cluster)[x_order],")"))+
  labs(x="",y="Hypoxia Score",title="SKCM_CRC mouse")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(color="black",size=10,angle = 30,hjust = 1,vjust = 1),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "top",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          strip.background = element_blank())
pdf('mouse_CD8Tcell_hypoxiaScore.pdf',5,3.5)
print(p)
dev.off()

###monocle
CD8Tcell@meta.data$functional.cluster <- as.character(CD8Tcell@meta.data$functional.cluster)
data_CD8Tcell <- as(as.matrix(CD8Tcell@assays$RNA@data), 'sparseMatrix')

pd_CD8Tcell <- new('AnnotatedDataFrame', data = CD8Tcell@meta.data)

fData_CD8Tcell <- data.frame(gene_short_name = row.names(data_CD8Tcell), row.names = row.names(data_CD8Tcell))
fd_CD8Tcell <- new('AnnotatedDataFrame', data = fData_CD8Tcell)

#Construct monocle cds
HSMM_CD8Tcell <- newCellDataSet(data_CD8Tcell,
                                       phenoData = pd_CD8Tcell,
                                       featureData = fd_CD8Tcell,
                                       lowerDetectionLimit = 0.5,
                                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150

#Run ordering algorithm
var_genes_CD8Tcell <- VariableFeatures(CD8Tcell)
ordering_genes_CD8Tcell <- var_genes_CD8Tcell

HSMM_CD8Tcell <- setOrderingFilter(HSMM_CD8Tcell, ordering_genes_CD8Tcell)
print(dim(exprs(HSMM_CD8Tcell)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM_CD8Tcell <- reduceDimension(HSMM_CD8Tcell,norm_method="none", 
                                        reduction_method="DDRTree",
                                        max_components=4,
                                        scaling=TRUE,
                                        verbose=TRUE,
                                        pseudo_expr=0)


## order cells change colors and theta to match your plot
HSMM_CD8Tcell <- orderCells(HSMM_CD8Tcell)
table(pData(HSMM_CD8Tcell)[,c('functional.cluster','State')])
HSMM_CD8Tcell <- orderCells(HSMM_CD8Tcell,root_state = 1) #set naive CD8T cell cluster as root

###color by pseudotime
pdf("mouse_CD8Tcell_Pseudotime.pdf",width = 4,height = 4)
plot_cell_trajectory(HSMM_CD8Tcell, 
                     theta = -10,
                     show_branch_points = F,
                     show_tree = TRUE, color_by = "Pseudotime", cell_size = 0.8) +
  scale_color_viridis_c()+
  labs(title = 'Pseudotime')+
  theme(legend.position = "bottom",legend.title = element_blank())
dev.off()


###color by hypoxia 
pdf("mouse_CD8Tcell_hypoxiaScore.pdf",width = 4,height = 4)
plot_cell_trajectory(HSMM_CD8Tcell, 
                     theta = -10,
                     show_branch_points = F,
                     show_tree = TRUE, color_by = "hypoxia", cell_size = 0.8) +
  scale_color_viridis_c()+
  labs(title = 'Hypoxia Score')+
  theme(legend.position = "bottom",legend.title = element_blank())
dev.off()

###color by cluster
pdf("3.figure/2.CD8Tcell_ClusterPseudotime.pdf",width = 4,height = 4)
plot_cell_trajectory(HSMM_CD8Tcell, 
                     color_by = "functional.cluster",
                     theta = -10,
                     show_branch_points = F,
                     show_tree = TRUE,
                     cell_size = 0.8) + 
  scale_color_manual(breaks=as.character(CD8TcellCluster$seurat_clusters), values=as.character(CD8TcellCluster$Color),name="") +
  theme(legend.position = "bottom")
dev.off()

