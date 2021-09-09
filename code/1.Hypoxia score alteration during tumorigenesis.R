###
setwd("/data/1.tumorigenesis")
options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(GSVA)

HypoxiaMarkers <- read.delim("HypoxiaMarkerSignature.txt",header=F)$V1
HypoxiaMarkersL <- list(HypoxiaMarkers=HypoxiaMarkers)

#load ID trans
GPL6480 <- read.delim("GPL6480_GSE17188.txt",comment.char = "#")[,c("ID","Gene.Symbol")]
colnames(GPL6480) <- c("ID","Gene.Symbol")
GPL6480 <- GPL6480[GPL6480$Gene.Symbol != "",]

#load expr matrix
GSE33479 <- read.delim("GSE33479_series_matrix.txt",comment.char = "!")
GSE33479 <- merge(GPL6480,GSE33479,by.x="ID",by.y = "ID_REF")
GSE33479 <- GSE33479[,-1]

#load Annotation
GSE33479_SampleAnn<- read.csv("GSE33479_SampleAnn.csv",header=F)  %>% dplyr::filter(V2 !="smoking status") %>%  tidyr::spread(V2,V3) 

StageClass <- data.frame(PhenotypicalStage = c("normal normofluorescent","normal hypofluorescent","hyperplasia","metaplasia","mild dysplasia","moderate dysplasia","severe dysplasia","carcinoma in situ","squamous cell carcinoma"),
                         Stage = paste("Stage",0:8,sep=""))
Stage0Sample <- GSE33479_SampleAnn[GSE33479_SampleAnn$`phenotypical stage`=="normal normofluorescent",]$V1
Stage8Sample <- GSE33479_SampleAnn[GSE33479_SampleAnn$`phenotypical stage`=="squamous cell carcinoma",]$V1


###Hypoxia expression across different stages
HypoxiaSigExp <- GSE33479[GSE33479$Gene.Symbol %in% HypoxiaMarkers,]
HypoxiaSigExp[,2:ncol(HypoxiaSigExp)] <- t(apply(HypoxiaSigExp[,2:ncol(HypoxiaSigExp)],1,scale))
HypoxiaSigExp <- HypoxiaSigExp[rev(order(HypoxiaSigExp$Gene.Symbol,apply(HypoxiaSigExp[,2:ncol(HypoxiaSigExp)],1,function(x){mean(x[names(x) %in%Stage8Sample])-mean(x[names(x) %in% Stage0Sample])}))),]
HypoxiaSigExp <- HypoxiaSigExp[!duplicated(HypoxiaSigExp$Gene.Symbol),]

HypoxiaSigExpM <- HypoxiaSigExp %>% tidyr::gather(key=GSM_ID,value=Exp,-Gene.Symbol) %>% tidyr::spread(Gene.Symbol,Exp) %>% dplyr::inner_join(GSE33479_SampleAnn[,c("phenotypical stage","V1")],by=c("GSM_ID"="V1")) %>%
  dplyr::inner_join(StageClass,by=c("phenotypical stage"="PhenotypicalStage"))
HypoxiaSigExpM$SumV <- apply(HypoxiaSigExpM[,2:16],1, sum)
HypoxiaSigExpM <- HypoxiaSigExpM[order(HypoxiaSigExpM$Stage,HypoxiaSigExpM$SumV),]
HypoxiaSigExp.m <- HypoxiaSigExp %>% tidyr::gather(key=GSM_ID,value=Exp,-Gene.Symbol) 
HypoxiaSigExp.m$Exp <- ifelse(HypoxiaSigExp.m$Exp >1.5,1.5,ifelse(HypoxiaSigExp.m$Exp < -1.5,-1.5,HypoxiaSigExp.m$Exp))
HypoxiaSigLast <- sapply(split(HypoxiaSigExpM,HypoxiaSigExpM$Stage),function(x){x[x$SumV==max(x$SumV),]$GSM_ID})

pdf("Hypoxia gene signature expression distribution.pdf",width = 6,height = 5)
ggplot(HypoxiaSigExp.m,aes(x=GSM_ID,y=Gene.Symbol,fill=Exp))+
  geom_tile(color=NA)+
  scale_x_discrete(limit=HypoxiaSigExpM$GSM_ID)+
  scale_fill_gradientn(limit=c(-1.5,1.5),colors= colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100))+
  theme(panel.background=element_rect(colour=NA,fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_text(color="black"),legend.background = element_blank(),
        axis.ticks = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
        strip.background = element_blank())+
  geom_vline(xintercept = match(HypoxiaSigLast,HypoxiaSigExpM$GSM_ID)+0.5,linetype="dashed")
dev.off()

###Hypoxia score in each stage
GSE33479.m <- as.matrix(GSE33479[grep("GSM",colnames(GSE33479),value=T)])
GSE33479.m <- t(apply(GSE33479.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE33479.m) <- GSE33479$Gene.Symbol
GSE33479_ES <- GSVA::gsva(GSE33479.m,HypoxiaMarkersL)

GSE33479_ES <- data.frame(t(GSE33479_ES))
GSE33479_ES$GSM_ID <-rownames(GSE33479_ES)

GSE33479_HypoImm <- GSE33479_ES %>% dplyr::inner_join(GSE33479_SampleAnn[,c("phenotypical stage","V1")],by=c("GSM_ID"="V1")) %>%
  dplyr::inner_join(StageClass,by=c("phenotypical stage"="PhenotypicalStage"))

GSE33479_HypoScore <- Rmisc::summarySE(GSE33479_HypoImm[,c("HypoxiaMarkers","Stage")], measurevar="HypoxiaMarkers", groupvars=c("Stage"))

pdf("Hypoxia score expression distribution.pdf",width = 6,height = 2.2)
ggplot(GSE33479_HypoScore,aes(x=Stage,y=HypoxiaMarkers,group=1))+
  geom_line(linetype = "dashed")+
  geom_errorbar(aes(ymin=HypoxiaMarkers-se,ymax=HypoxiaMarkers+se), width=0.2,position=position_dodge(0.9),size=0.5)+
  geom_point(color="orange",size=3) + 
  ylab("Hypoxia score")+
  scale_x_discrete(limit=GSE33479_HypoScore$Stage,labels=0:8)+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=10,color="black"),axis.title.y=element_text(size=10,color="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10,color="black"),axis.line = element_line(color="black"))
dev.off()  

