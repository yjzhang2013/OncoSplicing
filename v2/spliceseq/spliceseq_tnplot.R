library(ggplot2)
library(ggbeeswarm)
library(argparser, quietly=TRUE)

## Create a parser
p <- arg_parser("chunk size")
## Add command line arguments
p <- add_argument(p, "--ct", help="Cancer Type")
p <- add_argument(p, "--chunkN", help="chunk N")
p <- add_argument(p, "--chunkI", help="chunk I")
## Parse the command line arguments
argv <- parse_args(p)
## as.numeric transform
CancerType <- argv$ct
chunkN <- as.numeric(argv$chunkN)
chunkI <- as.numeric(argv$chunkI)

fileDir <- "/home/webdata/spliceseq/dataplot/"
outDir <- paste("/home/webdata/spliceseq/tnplots/",CancerType,"/",sep="")

if (file.exists(outDir)==FALSE){
  dir.create(outDir,recursive = T)
}

psi <- readRDS(paste(fileDir,"spliceseq_data_psi_",CancerType,".rds",sep=""))
aseAll <- names(psi)[names(psi)!="tissueType"]
aseAll <- aseAll[order(aseAll,decreasing = F)]
chunkSize <- length(aseAll)%/%chunkN
chunkLast <- chunkSize + length(aseAll)%%chunkN

start <- (chunkI-1)*chunkSize+1
if (chunkI<chunkN){
  end <- chunkI*chunkSize
}else if(chunkI==chunkN){
  end <- length(aseAll)
}

ase <- aseAll[start:end]

TNplot=function(SpliceEvent){
  cat("\n",SpliceEvent)
  outFilePath <- paste(outDir,CancerType,"-",SpliceEvent,"-TNplot.pdf",sep = "")
  
  df <- psi[,c("tissueType",SpliceEvent)]
  df$tissueType <- gsub("_","-",df$tissueType)
  df <- na.omit(df)
  tb <- data.frame(table(df$tissueType))
  names(df)[2] <- "PSI"
  if(nrow(tb)==2){
    width=3.5
    n <- paste(CancerType,"N",sep="-")
    t <- paste(CancerType,"T",sep="-")
    dfn <- df[df$tissueType==n,"PSI"]
    dft <- df[df$tissueType==t,"PSI"]
    ptest <- tryCatch(wilcox.test(dft,dfn)$p.value,error = function(e) return(NA))
    p <- format(ptest, scientific = TRUE,digits = 3)
    title=paste(CancerType,"-SpliceSeq\n",SpliceEvent,"\np-value = ",p,sep="")
    norm <- paste(n,"\nn=",length(dfn),sep="")
    tumor <- paste(t,"\nn=",length(dft),sep="")
    Xlimits <- c(n,t)
    Xlabels <- c(norm,tumor)
    color <- c("lightblue","#FF6A6A")
  }else{
    width=3
    t <- paste(CancerType,"T",sep="-")
    dft <- df[df$tissueType==t,"PSI"]
    title=paste(CancerType,"-SpliceSeq\n",SpliceEvent,sep="")
    tumor <- paste(t,"\nn=",length(dft),sep="")
    Xlimits <- c(t)
    Xlabels <- c(tumor)
    color <- c("#FF6A6A")
  }
  
  pdf(file=outFilePath,width=width,height=3.5)
  p <- ggplot(df,aes(x=tissueType,y=PSI,fill=tissueType)) +
    stat_boxplot(geom = "errorbar",width=0.15,color="grey20")+
    geom_boxplot(aes(x=tissueType,y=PSI),fill=color,outlier.fill="white",outlier.shape=21,outlier.size=2,width = 0.5)+
    scale_x_discrete(limits=Xlimits,labels = Xlabels)+
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
    labs(x="Tissue type",y="PSI value", title = title) +
    theme_bw()+
    theme(title=element_text(size=9),
          legend.position = "none",
          panel.grid.minor.y = element_line(color="white"),
          axis.text = element_text(size=9,color="black"),
          axis.title = element_text(size=9,face="bold"))
  print(p)
  dev.off()
}

sapply(ase,TNplot)
