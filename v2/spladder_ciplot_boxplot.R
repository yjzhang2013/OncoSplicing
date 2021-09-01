library(ggplot2)
library(ggbeeswarm)
library(argparser, quietly=TRUE)

## Create a parser
p <- arg_parser("chunk size")
## Add command line arguments
p <- add_argument(p, "--chunkN", help="chunk N")
p <- add_argument(p, "--chunkI", help="chunk I")
## Parse the command line arguments
argv <- parse_args(p)
## as.numeric transform
chunkN <- as.numeric(argv$chunkN)
chunkI <- as.numeric(argv$chunkI)


fileDir <- "/home/RNAseq/pancan/webdata/ClinAS/"
outDir <- "/home/RNAseq/pancan/webdata/ClinAS/plots/"
clin <- readRDS(file=paste(fileDir,"clinical_variates_with_sample_more_than20_spladder.rds",sep=""))
ls <- readRDS(paste(fileDir,"spladder_clinical_variates_as_PSIdata.rds",sep=""))
infoAS <- read.csv((paste(fileDir,"spladder_clinical_variates_as_info_10.csv",sep="")))
infoAS$infos <- paste(infoAS$Cancer_Type,infoAS$Gene_Symbol,infoAS$Splice_Event,infoAS$CI_Type,sep=";")

chunkSize <- nrow(infoAS)%/%chunkN
chunkLast <- chunkSize + nrow(infoAS)%%chunkN
infoAS$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))
infos <- infoAS[infoAS$chunk==chunkI,"infos"]

CIplot=function(infos){
  cat("\nchunkI-",chunkI,"-",infos)
  CancerType <- strsplit(infos,split=";")[[1]][1]
  GeneSymbol <- strsplit(infos,split=";")[[1]][2]
  SpliceEvent <- strsplit(infos,split=";")[[1]][3]
  CItype <- strsplit(infos,split=";")[[1]][4]
  
  outFilePath <- paste(outDir,CancerType,"-",GeneSymbol,"-",SpliceEvent,"-",CItype,".pdf",sep = "")
  ci <- clin[clin$cancerType==CancerType & clin$clType==CItype,]
  psi <- ls[[CancerType]]
  psi <- data.frame(t(psi))
  psi$sampleName <- row.names(psi)
  
  psi <- psi[,c("sampleName",SpliceEvent)]
  names(psi)[2] <- "PSI"
  df <- merge(ci,psi,by="sampleName")
  
  nPSI <- data.frame(table(df[!is.na(df$PSI),"clValue"]))
  names(nPSI) <- c("clValue","nValue")
  
  df1 <- df[df$clValue==nPSI$clValue[1],"PSI"]
  df2 <- df[df$clValue==nPSI$clValue[2],"PSI"]
  ptest <- tryCatch(wilcox.test(df1,df2)$p.value,error = function(e) return(NA))
  p <- format(ptest, scientific = TRUE,digits = 3)
  
  Xlimits <- nPSI$clValue
  Xlabels <- paste(nPSI$clValue,"\n(n=",nPSI$nValue,")",sep="")
  color <- c("lightblue","#FF6A6A")
  title <- paste("SplAdder - ",CancerType,"\n",GeneSymbol,"_",SpliceEvent,"\nP-value = ",p,sep="")
  
  pdf(file=outFilePath,width=3.5,height=3.5)
  p <- ggplot(df,aes(x=clValue,y=PSI,fill=clValue)) +
    stat_boxplot(geom = "errorbar",width=0.15,color="grey20")+
    geom_boxplot(aes(x=clValue,y=PSI),fill=color,outlier.fill="white",outlier.shape=21,outlier.size=2,width = 0.5)+
    scale_x_discrete(limits=Xlimits,labels = Xlabels)+
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
    labs(x=CItype,y="PSI value", title = title) +
    theme_bw()+
    theme(title=element_text(size=9),
          legend.position = "none",
          panel.grid.minor.y = element_line(color="white"),
          axis.text = element_text(size=9,color="black"),
          axis.title = element_text(size=9,face="bold"))
  suppressWarnings(print(p))
  dev.off()
}

sapply(infos,CIplot)
