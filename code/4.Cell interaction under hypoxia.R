setwd('/data/4.Cell interaction')
library(Seurat)
library(CellChat)
library(dplyr)
source("scores_RNA.r")

#load SKCM123139
GSE123139 <- readRDS('GSE123139_DefineTypes.rds')
Idents(GSE123139) <- 'Celltype..major.lineage.'
GSE123139 <- subset(GSE123139,idents = c('B','CD4Tconv','CD8Tex','DC','Mono/Macro','Plasma','Tprolif')) #subset immune cells

HypoxiaMarkers <- read.delim("../1.tumorigenesis/HypoxiaMarkerSignature.txt",header=F)$V1
HypoxiaMarkersL <- list(HypoxiaMarkers=HypoxiaMarkers)

HPscores <- score_cells(seur=GSE123139, names=HypoxiaMarkersL, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
GSE123139@meta.data$hypoxia <- HPscores[match(rownames(GSE123139@meta.data),rownames(HPscores)),'HypoxiaMarkers']

GSE123139M <- GSE123139@meta.data[grep("Mono",GSE123139@meta.data$Celltype..major.lineage.),]
GSE123139M$hypoxiaHL <- ifelse(GSE123139M$hypoxia > median(GSE123139M$hypoxia),'High','Low') #split by median
HighSample <- GSE123139M[grep('High',GSE123139M$hypoxiaHL),'Cell'] %>% as.character()
LowSample <- GSE123139M[grep('Low',GSE123139M$hypoxiaHL),'Cell'] %>% as.character()

GSE123139$Celltype..major.lineage. <- as.character(GSE123139$Celltype..major.lineage.)
GSE123139$Celltype..major.lineage.[HighSample] <- "Macro_HypoxiaHigh"
GSE123139$Celltype..major.lineage.[LowSample] <- "Macro_HypoxiaLow"

### create cellchat object 
cellchat <- createCellChat(object = GSE123139,group.by = 'Celltype..major.lineage.')

#LR data base 
CellChatDB <- CellChatDB.human

##pre-process
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat@LR$LRsig
cellchat <- projectData(cellchat,PPI.human)

##infer interaction
##liagnd receptor
cellchat <- computeCommunProb(cellchat,raw.use = FALSE,population.size = TRUE)
#filter cell<10
cellchat <- filterCommunication(cellchat,min.cells = 10)
df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name = 'netP')

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

#plot
levels(cellchat@idents)

vertex.receiver <- c(1:4,7,8)
pathway.show <- 'SPP1'

pdf('GSE123139_hierachyPlot_SPP1.pdf',8,4)
netVisual_aggregate(cellchat,signaling = pathway.show,vertex.receiver = vertex.receiver)
dev.off()


pdf('GSE123139_circle_SPP1_MIF.pdf',6,4)
par(mfrow = c(1,2), xpd=TRUE)
pathways.show <- c('SPP1')
pathways.show2 <- c('MIF')
netVisual_aggregate(cellchat,signaling = pathways.show,layout = 'circle')
netVisual_aggregate(cellchat,signaling = pathways.show2,layout = 'circle')
dev.off()
