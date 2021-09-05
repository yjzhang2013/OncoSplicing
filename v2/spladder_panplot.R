## plot SplAdder Pan-cancer TCGA project
## color modified in 2021.06.01
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
col1 <- rep(col1,16)[1:32]
col1[c(1,7,9,14,18,19,27,31,32)]<-"NA"
col2<- brewer.pal(8, "Dark2")[c(2,8)]
col2 <- rep(col2,16)[1:32]
colx <- c()
for(i in 1:32){
  colx <- c(colx,col1[i],col2[i])
}
colx<-colx[grepl("#",colx)]

colg <- brewer.pal(8, "Pastel2")[c(2,3)]
colg <- rep(colg,16)[1:31]
colors <- c(colx,colg)

## text
col1 <- brewer.pal(8, "Dark2")[c(2,8)]
col1 <- rep(col1,16)[1:32]
col1[c(1,7,9,14,18,19,27,31,32)]<-"NA"

colx <- c()
for(i in 1:32){
  colx <- c(colx,col1[i],col2[i])
}
colx<-colx[grepl("#",colx)]
colg <- rep("black",31)
colorXText <- c(colx,colg)

# def function
panplot <- function(ase){
  title=infoGDC[ase,"idType"]
  outFilePath <- paste("/home/RNAseq/pancan/panplot/",ase,".pdf",sep="")
  infoX <- infoAll[infoAll$SpliceEvent==ase,]
  
  dfp <- psiGDC[,c("tissueType",ase)]
  names(dfp)[2] <- "value"
  dfp[!gsub("-N","",gsub("-T","",dfp$tissueType))%in%infoX$tissueType,"value"] <- NA
  naSample <- row.names(dfp[is.na(dfp$value),])
  dfp$type <- "PSI value"
  
  dfi <- inReadGDC[,c("tissueType",ase)]
  dfi[naSample,ase] <- NA
  names(dfi)[2] <- "value"
  dfi$type <- "Reads-In"
  
  dfo <- outReadGDC[,c("tissueType",ase)]
  dfo[naSample,ase] <- NA
  names(dfo)[2] <- "value"
  dfo$type <- "Reads-Out"
  
  df <- rbind(dfp,dfi,dfo)
  df$type <- factor(df$type,levels = c("Reads-In","Reads-Out","PSI value"))
  
  ## the num of sample with value
  fclevel <- names(table(dfp$tissueType))
  nSample <- data.frame(table(dfp$tissueType))
  names(nSample) <- c("cancerType","nSample")
  nPSI <- data.frame(table(dfp[!is.na(dfp$value),"tissueType"]))
  names(nPSI) <- c("cancerType","nValue")
  nPSI <- merge(nSample,nPSI,by="cancerType",all.x=T)
  row.names(nPSI) <- nPSI$cancerType
  nPSI <- nPSI[fclevel,]
  # keep colors of tissue type with value
  colors <- colors[!is.na(nPSI$nValue) & nPSI$nValue!=0]
  
  nPSI[is.na(nPSI$nValue),"nValue"] <- 0
  tissueTypeN <- paste(nPSI$cancerType," (",nPSI$nValue,")",sep="")
  
  pdf(file=outFilePath,width=9,height=4.5,onefile = T)
  p <- ggplot(df,aes(x=tissueType,y=value,fill=tissueType)) +
    stat_boxplot(geom = "errorbar",width=0.3,color="grey20",size=0.2)+
    geom_boxplot(size=0.2,outlier.fill="white",outlier.shape=21,outlier.size=0.8,width = 0.55)+
    scale_x_discrete(limits=fclevel,breaks=fclevel,labels = tissueTypeN)+
    scale_fill_manual(values = colors)+
    labs(x="Tissue type",y="", title = title) +
    theme_bw()+
    facet_grid(type~.,scale="free",switch = "y") +
    geom_vline(xintercept = 55.5, size = 0.4, color="black",linetype = "dashed") +
    theme(title=element_text(size=8), 
          legend.position = "none",
          axis.text = element_text(size=6.7,color="black"),
          axis.title = element_text(size=8),
          axis.text.x = element_markdown(angle=90,hjust=1,vjust=0.5,color=colorXText),
          panel.grid = element_line(size=0.1),
          panel.grid.minor.y = element_blank(),
          #panel.grid.major.x = element_line(color=colp),
          panel.border = element_rect(size=0.2),
          strip.background = element_rect(color ="white",size = 0.5,fill = "white"), 
          strip.text = element_text(size=8),
          strip.placement = "outside")
  suppressWarnings(print(p))
  dev.off()
}

if (chunkI==0){
  ## prepare data
  infoGDC <- readRDS(file="/home/RNAseq/pancan/spladder_infoAllTis.rds")
  infoGDC <- infoGDC[order(infoGDC$idType),]
  names(infoGDC)
  colName <- c("SpliceEvent","GeneSymbol","SpliceType","idType")
  infoGDC <- infoGDC[,colName]
  aseAll <- infoGDC[!duplicated(infoGDC$SpliceEvent),c("SpliceEvent","idType")]
  aseAll$order <- seq(1:nrow(aseAll))
  row.names(aseAll) <- aseAll$SpliceEvent
  
  chunkSize <- nrow(aseAll)%/%chunkN
  chunkLast <- chunkSize + nrow(aseAll)%%chunkN
  aseAll$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))
  gr.info <- aseAll
  gr.info <- split(gr.info,gr.info$chunk)
  gr.psi <- readRDS(file="/home/RNAseq/pancan/dataPan/tcga_gtex_psi_tissueType.rds")
  gr.inRead <- readRDS(file="/home/RNAseq/pancan/dataPan/tcga_gtex_inRead_tissueType.rds")
  gr.outRead <- readRDS(file="/home/RNAseq/pancan/dataPan/tcga_gtex_outRead_tissueType.rds")
  
  for (i in 1:chunkN){
    infoGDC <- gr.info[[i]]
    psiGDC <- gr.psi[,c("tissueType",row.names(infoGDC))]
    inReadGDC <- gr.inRead[,c("tissueType",row.names(infoGDC))]
    outReadGDC <- gr.outRead[,c("tissueType",row.names(infoGDC))]
    save(infoGDC,psiGDC,inReadGDC,outReadGDC,file=paste("/home/RNAseq/pancan/gdcAS/tcga3/temp_psi_read_info_fix_",i,".Rdata",sep=""),compress = T)
  }
}else{
  ## load data and plot
  load(file=paste("/home/RNAseq/pancan/temp_psi_read_info_fix_",chunkI,".Rdata",sep=""))
  infoAll <- readRDS(file="/home/RNAseq/pancan/spladder_infoAllTis_tcga_gtex.rds")
  ase <- infoGDC$SpliceEvent
  #aseDone <- readRDS(file="/home/RNAseq/pancan/asePanPlotDone.rds")
  #ase <- ase[!ase%in%aseDone]
  sapply(ase,panplot)
}
