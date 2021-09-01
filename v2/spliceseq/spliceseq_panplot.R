## plot SpliceSeq Pan-cancer TCGA project
library(RColorBrewer)
library(ggplot2)
library(ggtext)
library(argparser, quietly=TRUE)

## Create a parser
p <- arg_parser("chunk size")
## Add command line arguments
p <- add_argument(p, "--chunkN", help="chunk number")
p <- add_argument(p, "--chunkI", help="chunk index")
## Parse the command line arguments
argv <- parse_args(p)
## as.numeric transform
chunkN <- as.numeric(argv$chunkN)
chunkI <- as.numeric(argv$chunkI)


## color management
col1 <- brewer.pal(8, "Pastel2")[c(2,8)]
col1 <- rep(col1,17)[1:33]
col1[c(1,7,14,15,19,20,28,32,33)]<-"NA"
col2<- brewer.pal(8, "Dark2")[c(2,8)]
col2 <- rep(col2,17)[1:33]
colx <- c()
for(i in 1:33){
  colx <- c(colx,col1[i],col2[i])
}
colx<-colx[grepl("#",colx)]
colors <- colx


col1 <- brewer.pal(8, "Dark2")[c(2,8)]
col1 <- rep(col1,17)[1:33]
col1[c(1,7,14,15,19,20,28,32,33)]<-"NA"
colx <- c()
for(i in 1:33){
  colx <- c(colx,col1[i],col2[i])
}
colx<-colx[grepl("#",colx)]
colorXText <- colx


panplot <- function(ase){
  dfp <- psiGDC[,c("tissueType",ase)]
  names(dfp)[2] <- "value"
  
  ## the #sample with value, paste as XText 
  fclevel <- names(table(dfp$tissueType))
  nSample <- data.frame(table(dfp$tissueType))
  names(nSample) <- c("cancerType","nSample")
  nPSI <- data.frame(table(dfp[!is.na(dfp$value),"tissueType"]))
  names(nPSI) <- c("cancerType","nValue")
  nPSI <- merge(nSample,nPSI,by="cancerType",all.x=T)
  row.names(nPSI) <- nPSI$cancerType
  nPSI <- nPSI[fclevel,]
  # keep color of tissue type with value
  colors <- colors[!is.na(nPSI$nValue)]
  
  nPSI[is.na(nPSI$nValue),"nValue"] <- 0
  #tissueTypeN <- paste(nPSI$cancerType," (",nPSI$nValue,"/",nPSI$nSample,")",sep="")
  tissueTypeN <- paste(nPSI$cancerType," (",nPSI$nValue,")",sep="")
  
  title=ase
  outFilePath <- paste("/home/u1357/RNAseq/pancan/spliceSeq/panplot/fix/",ase,".pdf",sep="")
  pdf(file=outFilePath,width=9,height=3.5,onefile = T)
  p<-ggplot(dfp,aes(x=tissueType,y=value,fill=tissueType)) +
    stat_boxplot(geom = "errorbar",width=0.3,color="grey20",size=0.2)+
    geom_boxplot(size=0.2,outlier.fill="white",outlier.shape=21,outlier.size=1,width = 0.55)+
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
    scale_x_discrete(limits=fclevel,breaks=fclevel,labels = tissueTypeN)+
    scale_fill_manual(values = colors)+
    labs(x="Tissue type",y="PSI value", title = title) +
    theme_bw()+
    theme(title=element_text(size=9),
          legend.position = "none",
          axis.text = element_text(size=8.5,color="black"),
          axis.title = element_text(size=9),
          axis.text.x = element_markdown(angle=90,hjust=1,vjust=0.5,color=colorXText),
          panel.grid = element_line(size=0.1)
    )
  suppressWarnings(print(p))
  dev.off()
}


if (chunkI==0){
  ## prepare data
  infoGDC <- readRDS(file="/home/u1357/RNAseq/pancan/spliceSeq/infoNew/infoAllTis.rds")
  #aseAll <- infoGDC %>% group_by(SpliceEvent,idType) %>% summarise(n=n())
  aseAll <- infoGDC[!duplicated(infoGDC$SpliceEvent),c("SpliceEvent","idType")]
  aseAll <- aseAll[order(aseAll$idType),]
  aseAll$order <- seq(1:nrow(aseAll))
  
  chunkSize <- nrow(aseAll)%/%chunkN
  chunkLast <- chunkSize + nrow(aseAll)%%chunkN
  aseAll$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))
  saveRDS(aseAll,file="/home/u1357/RNAseq/pancan/spliceSeq/psiData/aseAll_panplot_psi.rds")
  
  gr.info <- readRDS(file="/home/u1357/RNAseq/pancan/spliceSeq/psiData/aseAll_panplot_psi.rds")
  gr.psi <- readRDS(file="./data_psi_33_cancer_types.rds")
  gr.info <- split(gr.info,gr.info$chunk)
  #i=1
  for (i in 1:chunkN){
    infoGDC <- gr.info[[i]]
    psiGDC <- gr.psi[,c("tissueType",infoGDC$SpliceEvent)]
    save(infoGDC,psiGDC,file=paste("temp_psi_info_",i,".Rdata",sep=""),compress = T)
  }
}else{
  ## load data and plot
  load(file=paste("temp_psi_info_",chunkI,".Rdata",sep=""))
  ase <- infoGDC$SpliceEvent
  sapply(ase,panplot)
}
