setwd('/data/6.ICB')
library(dplyr)
library(ggplot2)
library(ggpubr)
HypoxiaMarkers <- read.delim("../1.tumorigenesis/HypoxiaMarkerSignature.txt",header=F)$V1
HypoxiaMarkersL <- list(HypoxiaScore=HypoxiaMarkers)

##load
IMvigor210 <- readRDS('IMvigor210.rds')

IMvigor210m <- as.matrix(IMvigor210$gset_data_matrix_symbol[[1]])
IMvigor210_gsva <- GSVA::gsva(IMvigor210m,HypoxiaMarkersL,method='gsva')
IMvigor210_gsva <- t(IMvigor210_gsva) %>% as.data.frame()
IMvigor210_gsva$sample_use <- rownames(IMvigor210_gsva)

IMvigor210_meta <- IMvigor210$gset_data_pData[[1]]
IMvigor210_all <- merge(IMvigor210_gsva,IMvigor210_meta,by = 'sample_use')
IMvigor210_allF <- IMvigor210_all[-grep('None',IMvigor210_all$response_use),]

compare <- list(c("Benefit","NonBenefit"))
p <- ggplot(IMvigor210_allF,aes(x=response_use,y=HypoxiaScore))+
  geom_boxplot(aes(fill=response_use),width = 0.6)+
  scale_fill_manual(values = c('#456AAA','#C21B21'))+
  stat_compare_means(comparisons = compare,method = "t.test") + #,label = "p.signif"
  labs(x="",y="Hypoxia Score",title="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "top",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_blank(),
          strip.background = element_blank())
print(p)
