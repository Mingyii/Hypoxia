setwd("/data/3.ImmuneAI")
my.cor.test <- function(...) {
  obj<-try(cor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
library(dplyr)
library(magrittr)
library(ggplot2)

#load TCGA expr
HypoxiaScore.Exp <- readr::read_rds("pancan33_HypoxiaScore.rds")

ImmuCellAI <- readr::read_rds("TCGA_ImmuCellAI.rds.gz")
ImmuCellAI <- ImmuCellAI %>% dplyr::mutate(AI=purrr::map(.x=AI,function(.x){
  .x %>% dplyr::mutate(barcode=rownames(.))
}))

HypoxiaImmuCellAIScore <- HypoxiaScore.Exp %>% dplyr::inner_join(ImmuCellAI,by="cancer_types") %>%
  dplyr::mutate(HypoxiaImmSig = purrr::map2(.x = HypoxiaScore,.y=AI,function(.x,.y){
    dplyr::inner_join(.x,.y,by="barcode")
  }))

HypoxiaImmuCellAIScore <- HypoxiaImmuCellAIScore %>% dplyr::mutate(Corr = purrr::map(.x = HypoxiaImmSig ,function(.x){
  .x <- .x %>% dplyr::filter(substr(barcode,14,15) %in% c("01","03")) 
  Corr = data.frame(ImmuneCells = colnames(.x)[3:ncol(.x)] )
  Corr[,c("estimate.cor","p.value")] <- t(apply(.x[,as.character(Corr$ImmuneCells)],2,function(x)unlist(my.cor.test(as.numeric(x),as.numeric(.x$HypoxiaScore),method="spearman"))))
  Corr$FDR <- p.adjust(Corr$p.value,method="fdr")
  Corr$Class <- ifelse(Corr$estimate.cor > 0 & Corr$FDR < 0.05,1,ifelse(Corr$estimate.cor < 0 & Corr$FDR < 0.05,-1,0))
  return(Corr)
})) 


HypoxiaImmuCellAIScoreAll <- dplyr::bind_rows(HypoxiaImmuCellAIScore$Corr) %>%
  dplyr::mutate(cancer_types = rep(HypoxiaImmuCellAIScore$cancer_types,times=unlist(lapply(HypoxiaImmuCellAIScore$Corr,nrow))) )

HypoxiaImmuCellAIScoreCount <- t(sapply(split(HypoxiaImmuCellAIScoreAll[,"Class"],HypoxiaImmuCellAIScoreAll$ImmuneCells),function(x)table(factor(x,levels=c(-1,0,1))))) %>%
  data.frame() %>% set_colnames(c("NegCor","NonCor","PosCor")) %>% dplyr::mutate(ImmuneCells=row.names(.))
HypoxiaImmuCellAIScoreCount <- HypoxiaImmuCellAIScoreCount[order(HypoxiaImmuCellAIScoreCount$PosCor - HypoxiaImmuCellAIScoreCount$NegCor),]
HypoxiaImmuCellAIScoreCancers <- t(sapply(split(HypoxiaImmuCellAIScoreAll[,"Class"],HypoxiaImmuCellAIScoreAll$cancer_types),function(x)table(factor(x,levels=c(-1,0,1))))) %>%
  data.frame() %>% set_colnames(c("NegCor","NonCor","PosCor")) %>% dplyr::mutate(cancer_types=rownames(.))
HypoxiaImmuCellAIScoreCancers <- HypoxiaImmuCellAIScoreCancers[order(HypoxiaImmuCellAIScoreCancers$PosCor - HypoxiaImmuCellAIScoreCancers$NegCor),]
HypoxiaImmuCellAIScoreAllm <- HypoxiaImmuCellAIScoreAll[HypoxiaImmuCellAIScoreAll$ImmuneCells !="InfiltrationScore",]
y_lab <- HypoxiaImmuCellAIScoreCount$ImmuneCells[HypoxiaImmuCellAIScoreCount$ImmuneCells !="InfiltrationScore"]

pdf("HypoxiaScore_ImmuCellAI_Cor.pdf",width=7,height = 5)
ggplot(HypoxiaImmuCellAIScoreAllm)+
  geom_point(aes(x =cancer_types,y = ImmuneCells,color = estimate.cor,size=-log10(FDR+10^-10)))+
  scale_x_discrete(limit=HypoxiaImmuCellAIScoreCancers$cancer_types)+
  scale_size_continuous(limits = c(0,10),range=c(1,4),breaks=c(2,5,10),labels=c(0.01,10^-5,10^-10),name="FDR")+
  scale_color_gradientn(limit=c(-0.8,0.8),colors= colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100),breaks=seq(-0.5,0.5,length.out = 3),labels=seq(-0.5,0.5,length.out = 3),name="Spearman Corr.")+
  scale_y_discrete(limit=y_lab,labels=gsub("\\.","\\ ",y_lab))+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title=element_blank(),axis.text.x=element_text(color="black",angle=90,vjust=0.5,hjust=1),
          axis.text.y=element_text(color="black"),legend.background = element_blank(),
          axis.ticks = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          strip.background = element_blank())+
  geom_point(data = HypoxiaImmuCellAIScoreAllm[HypoxiaImmuCellAIScoreAllm$FDR < 0.05,],aes(x =cancer_types,y = ImmuneCells,size=-log10(FDR+10^-10)),shape=21,color="black")
dev.off()


HypoxiaImmuCellAIScoreCountm <- HypoxiaImmuCellAIScoreCount %>% tidyr::gather(key=Cor,value=Count,-ImmuneCells)

pdf("HypoxiaScore_ImmuCellAI_CancerNumbers.pdf",width=3,height = 5)
ggplot(HypoxiaImmuCellAIScoreCountm[HypoxiaImmuCellAIScoreCountm$Cor %in% c("NegCor","PosCor"),],aes(x=ImmuneCells,fill=Cor,y=Count))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_x_discrete(limit=y_lab,labels=gsub("\\.","\\ ",y_lab),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(limits=c("NegCor","PosCor"),values=c("blue","red"),labels=c("Neg","Pos"))+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title=element_blank(),axis.text.x=element_text(color="black"),
          axis.text.y=element_text(color="black"),legend.background = element_blank(),axis.ticks.y=element_blank(),
          axis.line.x = element_line(color = "black"),legend.position = c(0.8,0.9))
dev.off()