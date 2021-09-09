setwd('/data/2.spatial')
library(Seurat)
library(ggpubr)
library(ggplot2)
library(dplyr)
source("scores_spatial.r")

#CRC patient1 spatial
###Cluster and hypxoia cluster boxplot

#load data
CRC_ST1 <- readRDS('CRC_ST1_DefineTypes.rds')
Idents(CRC_ST1) <- 'DefineTypes'

pdf('spatial/3.figure/1.DimPlot_ST1.pdf',5,7)
SpatialDimPlot(CRC_ST1) +theme(legend.position = 'bottom')
dev.off()

HypoxiaMarkers <- read.delim("../1.tumorigenesis/HypoxiaMarkerSignature.txt",header=F)$V1
HypoxiaMarkersL <- list(HypoxiaMarkers=HypoxiaMarkers)

HPscores <- score_cells(seur=CRC_ST1, names=HypoxiaMarkersL, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
CRC_ST1@meta.data$Hyopxia_score <- HPscores[match(rownames(CRC_ST1@meta.data),rownames(HPscores)),]
CRC_ST1_meta <- CRC_ST1@meta.data

x_order <- lapply(split(CRC_ST1_meta$Hyopxia_score,CRC_ST1_meta$DefineTypes), median) %>% unlist()
x_order <- names(x_order)[order(x_order)]

##box plot
pdf("BoxPlot_CRC_ST1_hypoxiaScore.pdf",width=6,height=3.5)
ggplot(CRC_ST1_meta,aes(x=DefineTypes,y=Hyopxia_score))+
  geom_boxplot(aes(fill=DefineTypes),width = 0.6)+
  scale_x_discrete(limits=x_order,label=paste0(x_order,"(",table(CRC_ST1_meta$DefineTypes)[x_order],")"))+
  labs(x="",y="Hypoxia Score",title="CRC/ST1")+
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
dev.off()