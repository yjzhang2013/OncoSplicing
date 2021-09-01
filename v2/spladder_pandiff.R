library(gridExtra)
library(ggrepel)
library(ggplot2)
library(scales)
library(ggtext)
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

fileDir <- "/home/RNAseq/pancan/webdata/PanCancer/SplAdder/panDataPlot/data/"
outDir <- "/home/RNAseq/pancan/webdata/PanCancer/SplAdder/panDataPlot/plots/diff/"

infoAS <- readRDS(file=paste(fileDir,"spladder_infoAllTis_all.rds",sep=""))

aseAll <- infoAS[!duplicated(infoAS$SpliceEvent),c("SpliceEvent","idType")]
aseAll$order <- seq(1:nrow(aseAll))
chunkSize <- nrow(aseAll)%/%chunkN
chunkLast <- chunkSize + nrow(aseAll)%%chunkN
aseAll$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))

ase <- aseAll[aseAll$chunk==chunkI,"SpliceEvent"]

infoAS <- infoAS[infoAS$SpliceEvent%in%ase,]


# def PanDiff plot
PanDiff <- function(SpliceEvent){
  infoX <- infoAS[infoAS$SpliceEvent==SpliceEvent,]
  idType <- aseAll[aseAll$SpliceEvent==SpliceEvent,"idType"]
  outFilePath <- paste(outDir,SpliceEvent,"-PanDiff.pdf",sep = "")
  
  ## diff
  dfncol <- c("SpliceEvent","cancerType","NdtPSI","medtPSI","NdnPSI","mednPSI","psiDiff","bhDiff","idType")
  dfgcol <- c("SpliceEvent","cancerType","NdtPSI","medtPSI","NdgPSI","medgPSI","psiDiffGTEx","bhDiffGTEx","idType")
  
  dfn <- infoX[,dfncol]
  dfn <- na.omit(dfn)
  dfn$logbhDiff <- -log10(dfn$bhDiff)
  
  dfg <- infoX[,dfgcol]
  names(dfg) <- dfncol
  dfg <- na.omit(dfg)
  dfg$logbhDiff <- -log10(dfg$bhDiff)
  
  difplot <- function(df,type,colors){
    title <- idType
    if (type=="TN"){
      Xtitle <- "PSI difference (Tumor-Normal)"
    }else if(type=="TG"){
      Xtitle <- "PSI difference (Tumor-GTEx)"
    }
    ggplot(df, aes(psiDiff,logbhDiff,fill=cancerType,label = cancerType)) +
      geom_text_repel(size=3) +
      geom_point(aes(size = medtPSI),shape = 21,  stroke=0.2, alpha=0.5, show.legend = T) +
      scale_size_continuous(limits=sizelimits, breaks=sizebreaks)+
      theme_bw(base_size = 9) +
      guides(fill = FALSE,size=guide_legend(title="Tumor PSI")) +
      scale_fill_manual(values=colors) +
      scale_y_continuous(limits=Ylimits,breaks=Ybreaks,labels=Ylabels,position="left") +
      scale_x_continuous(limits=Xlimits,breaks=Xbreaks,labels=Xlabels,position="bottom") +
      labs(x=Xtitle,y="-log10(FDR)", title = title) +
      geom_hline(yintercept = cutoff, size = 0.4, color="#dd7d6a",linetype = "dashed") +
      theme(title=element_text(size=9),
            #axis.text.x = element_text(size=9,color=Xcolor),
            #axis.text.y = element_text(size=9,color=Ycolor),
            axis.text.x = element_markdown(size=9,color=Xcolor),
            axis.text.y = element_markdown(size=9,color=Ycolor),
            axis.title = element_text(size=9),
            legend.title = element_text(size=9))
  }
  
  ## X-axis
  df <- rbind(dfn,dfg)
  maxdif <- max(abs(df$psiDiff))
  n = floor(-log(maxdif, 10) + 1)
  maxval <- round(maxdif,n)
  Xlimits <- c(-max(maxdif,maxval),max(maxdif,maxval))
  Xbreaks <- c(-maxval,-maxval/2,0,maxval/2,maxval)
  Xquat <- maxval/2
  
  ## def function for x-axis transform
  Xtransf <- function(data){
    xCut <- 0.2
    if(abs(data)<=xCut){
      re <- data*Xquat/xCut
    }else if(data>xCut){
      re <- (1+(data-xCut)/(maxval-xCut))*Xquat
    }else if(data< -xCut){
      re <- (-1+(data-(-xCut))/(maxval-xCut))*Xquat
    }
    return(re)
  }
  
  Xtf <- (length(dfn$psiDiff[dfn$psiDiff<=0.2])>1 | length(dfg$psiDiff[dfg$psiDiff<=0.2])>1)
  
  if (maxval>0.8 & Xtf){
    dfn$psiDiff <- sapply(dfn$psiDiff,Xtransf)
    dfg$psiDiff <- sapply(dfg$psiDiff,Xtransf)
    Xlabels <- c(-maxval,-0.2,0,0.2,maxval)
    Xcolor <- c("black","red","black","red","black")
  }else{
    Xlabels <- Xbreaks
    Xcolor <- rep("black",length(Xbreaks))
  }
  
  
  ## Y-axis
  df$logbhDiff <- -log10(df$bhDiff)
  maxbh <- max(df$logbhDiff)
  y = floor(-log(maxbh, 10) + 1)
  #Ymaxval <- round(maxbh,y)
  Ymaxval <- ifelse(round(maxbh,y)<2,2,round(maxbh,y)) ## 2021.05.31
  Ylimits <- c(0,max(maxbh,Ymaxval))
  Ybreaks <- c(0,Ymaxval*0.25,Ymaxval*0.5,Ymaxval*0.75,Ymaxval)
  Yquat <- Ymaxval*0.25
  
  Ytransf <- function(data){
    yCut <- 2
    if(data<=yCut){
      re <- data*Yquat/yCut
    }else if(data>yCut){
      re <- 1*Yquat+((data-yCut)/(Ymaxval-yCut))*3*Yquat
    }
    return(re)
  }
  
  Ytf <- (length(dfn$logbhDiff[dfn$logbhDiff<=2])>1 | length(dfg$logbhDiff[dfg$logbhDiff<=2])>1)
  
  if (Ymaxval>20 & Ytf){
    dfn$logbhDiff <- sapply(dfn$logbhDiff,Ytransf)
    dfg$logbhDiff <- sapply(dfg$logbhDiff,Ytransf)
    Ylabels <- c(0,2,round((Ymaxval-2)/3,1)+2,round(2*(Ymaxval-2)/3,1)+2,Ymaxval)
    Ycolor <- c("black","red","black","black","black")
    cutoff <- sapply(-log(0.05,10),Ytransf)
  }else{
    Ylabels <- Ybreaks
    Ycolor <- rep("black",length(Ybreaks))
    cutoff <- -log(0.05,10)
  }
  
  ## point size - tumor psi
  maxsize <- max(df$medtPSI)
  minsize <- min(df$medtPSI)
  range <- maxsize-minsize
  
  if (range==0){
    sizelimits = c(minsize,maxsize)
    sizebreaks = c(minsize,maxsize)
  }else if (range<0.05){
    sizelimits = c(minsize,maxsize)
    sizebreaks = c(minsize,maxsize)
  }else{
    x = ceiling(-log(maxsize, 10) + 1)
    n = ceiling(-log(minsize, 10) + 1)
    maxval <- round(maxsize,x)
    minval <- round(minsize,n)
    sizelimits = c(min(minval,minsize),max(maxval,maxsize))
    if(minval==0){
      sizebreaks <- round(seq(minval,maxval,(maxval-minval)/4)[-1],2)
    }else{
      sizebreaks <- round(seq(minval,maxval,(maxval-minval)/4),2)
    }
  }
  
  ## color
  dfcol <- data.frame(table(df$cancerType))
  if(nrow(dfcol)>0){
    colors <-  hue_pal()(nrow(dfcol))
    colors1 <- colors[dfcol$Var1%in%dfn$cancerType]
    colors2 <- colors[dfcol$Var1%in%dfg$cancerType]
  }
  
  nCancer=3
  if (nrow(dfn)>=nCancer & nrow(dfg)>=nCancer){
    width=10
    p1 <- difplot(dfn,"TN",colors1)
    p2 <- difplot(dfg,"TG",colors2)
    
    pdf(file=outFilePath,width=width,height=4)
    suppressWarnings(grid.arrange(p1,p2, ncol=2))
    dev.off()
  }else if (nrow(dfn)>=nCancer & nrow(dfg)<nCancer){
    width=5
    p1 <- difplot(dfn,"TN",colors1)
    pdf(file=outFilePath,width=width,height=4)
    suppressWarnings(print(p1))
    dev.off()
  }else if (nrow(dfn)<nCancer & nrow(dfg)>=nCancer){
    width=5
    p2 <- difplot(dfg,"TG",colors2)
    pdf(file=outFilePath,width=width,height=4)
    suppressWarnings(print(p2))
    dev.off()
  }else if (nrow(dfn)<nCancer & nrow(dfg)<nCancer){
    library(gplots)
    width=5
    pdf(file=outFilePath,width=width,height=4)
    temptext1 <- paste("! ",SpliceEvent,"\nNo data presented.\n","Try another one!",sep="")
    textplot(temptext1,valign="top", cex=1.3, halign= "center",col="grey")
    dev.off()
  }
  return(outFilePath)
}

## plot
sapply(ase,PanDiff)
