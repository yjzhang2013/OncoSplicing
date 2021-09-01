library(gridExtra)
library(ggrepel)
library(ggplot2)
library(scales)
library(ggtext)
library(gplots)
library(argparser, quietly=TRUE)

## Create a parser
p <- arg_parser("chunk size")
## Add command line arguments
p <- add_argument(p, "--st", help="survival type")
p <- add_argument(p, "--chunkN", help="chunk N")
p <- add_argument(p, "--chunkI", help="chunk I")
## Parse the command line arguments
argv <- parse_args(p)
## as.numeric transform
survType <- argv$st
chunkN <- as.numeric(argv$chunkN)
chunkI <- as.numeric(argv$chunkI)

fileDir <- "/home/webdata/spladder/"
outDir <- "/home/webdata/spladder/pancox/"

# file include each survival data of all cancer types
fileName <- paste(fileDir,"spladder_ase_all_cancertype_",survType,".rds",sep="")
infoAS <- readRDS(file=fileName)

aseAll <- read.csv(file="/home/u1357/webdata/spladder/newAdd/spladder_pancan_info.csv")
aseAll$idType <- paste(aseAll$Gene_Symbol,aseAll$Splice_Event,sep="_")
aseAll <- aseAll[,c("Splice_Event","idType")]
names(aseAll) <- c("SpliceEvent","idType")
aseAll <- aseAll[order(aseAll$SpliceEvent,decreasing = F),]
chunkSize <- nrow(aseAll)%/%chunkN
chunkLast <- chunkSize + nrow(aseAll)%%chunkN
aseAll$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))

aseAll <- aseAll[aseAll$chunk==chunkI,]
ase <- aseAll$SpliceEvent
infoAS <- infoAS[infoAS$SpliceEvent%in%ase,]
panx <- paste("Pan",survType,sep="")

# def function
PanCox <- function(SpliceEvent){
  cat("\n",SpliceEvent)
  infoX <- infoAS[infoAS$SpliceEvent==SpliceEvent,]
  idType <- aseAll[aseAll$SpliceEvent==SpliceEvent,"idType"]
  outFilePath <- paste(outDir,SpliceEvent,"-",panx,".pdf",sep = "")
  
  ## coxph
  dfmcol <- c("SpliceEvent","cancerType","cutmed","pvalHRmed","bhHRmed","HRmed","nMinHRmed","nEventHRmed","idType")
  dffcol <- c("SpliceEvent","cancerType","cutfit","pvalHRfit","bhHRfit","HRfit","nMinHRfit","nEventHRfit","idType")
  
  dfm <- infoX[,dfmcol]
  dfm <- na.omit(dfm)
  dfm <- dfm[dfm$nMinHRmed>10 & dfm$nEventHRmed>5,]
  dfm$logHRmed <- log2(dfm$HRmed)
  dfm$logpvalHRmed <- ifelse(dfm$pvalHRmed==0,17,-log10(dfm$pvalHRmed))
  
  dff <- infoX[,dffcol]
  names(dff) <- dfmcol
  dff <- na.omit(dff)
  dff <- dff[dff$nMinHRmed>10 & dff$nEventHRmed>5,]
  dff$logHRmed <- log2(dff$HRmed)
  dff$logpvalHRmed <- ifelse(dff$pvalHRmed==0,17,-log10(dff$pvalHRmed))
  
  ## def plot
  coxplot <- function(df,type,colors){
    title <- idType
    if (type=="med"){
      Xtitle <- paste("log2(Hazard Ratio), ",survType,"\nMedian cutoff",sep="")
    }else if(type=="fit"){
      Xtitle <- paste("log2(Hazard Ratio), ",survType,"\nOptimal cutoff",sep="")
    }
    ggplot(df, aes(logHRmed,logpvalHRmed,fill=cancerType,label = cancerType)) +
      geom_text_repel(size=3) +
      geom_point(aes(size = cutmed),shape = 21,  stroke=0.2, alpha=0.5, show.legend = T) +
      scale_size_continuous(limits=sizelimits, breaks=sizebreaks)+
      theme_bw(base_size = 9) +
      guides(fill = F,size=guide_legend(title="Cutoff")) +
      scale_fill_manual(values=colors) +
      scale_y_continuous(limits=Ylimits,breaks=Ybreaks,labels=Ylabels,position="left") +
      scale_x_continuous(limits=Xlimits,breaks=Xbreaks,labels=Xlabels,position="bottom") +
      labs(x=Xtitle,y="-log10(P-value)", title = title) +
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
  df <- rbind(dfm,dff)
  df$logHRmed <- log2(df$HRmed)
  maxhr <- max(abs(df$logHRmed))
  n = floor(-log(maxhr, 10) + 1)
  maxval <- round(maxhr,n)
  Xlimits <- c(-max(maxhr,maxval),max(maxhr,maxval))
  Xbreaks <- c(-maxval,-maxval/2,0,maxval/2,maxval)
  Xquat <- maxval/2
  
  # def function for x-axis transform
  Xtransf <- function(data){
    xCut <- 2
    if(abs(data)<=xCut){
      re <- data*Xquat/xCut
    }else if(data>xCut){
      re <- (1+(data-xCut)/(maxval-xCut))*Xquat
    }else if(data< -xCut){
      re <- (-1+(data-(-xCut))/(maxval-xCut))*Xquat
    }
    return(re)
  }
  
  Xtf <- length(dfm$logHRmed[dfm$logHRmed<=2])>1 | length(dff$logHRmed[dff$logHRmed<=2])>1
  
  if (maxval>10 & Xtf){
    dfm$logHRmed <- sapply(dfm$logHRmed,Xtransf)
    dff$logHRmed <- sapply(dff$logHRmed,Xtransf)
    Xlabels <- c(-maxval,-2,0,2,maxval)
    Xcolor <- c("black","red","black","red","black")
  }else{
    Xlabels <- Xbreaks
    Xcolor <- rep("black",length(Xbreaks))
  }
    
    
  ## Y-axis
  # min(p)<10-16
  df$logpvalHRmed <- ifelse(df$pvalHRmed==0,16,-log10(df$pvalHRmed))
  maxbh <- max(df$logpvalHRmed)
  y = floor(-log(maxbh, 10) + 1)
  Ymaxval <- ifelse(round(maxbh,y)<2,2,round(maxbh,y))
  Ylimits <- c(0,max(maxbh,Ymaxval))
  
  if (maxbh==16){
    Ybreaks <- c(0,Ymaxval*0.25,Ymaxval*0.5,Ymaxval*0.75,16,Ymaxval)
    Ylabels <- c(Ybreaks[1:5],"Inf")
    Ycolor <- c("black","black","black","black","black","red")
    cutoff <- -log(0.05,10)
  }else{
    Ybreaks <- c(0,Ymaxval*0.25,Ymaxval*0.5,Ymaxval*0.75,Ymaxval)
    Ylabels <- Ybreaks
    Ycolor <- rep("black",length(Ybreaks))
    cutoff <- -log(0.05,10)
  }
  
  
  ## point size - psi cutoff
  maxsize <- max(df$cutmed)
  minsize <- min(df$cutmed)
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
    colors1 <- colors[dfcol$Var1%in%dfm$cancerType]
    colors2 <- colors[dfcol$Var1%in%dff$cancerType]
  }
  
  nCancer=3
  if (nrow(dfm)>=nCancer & nrow(dff)>=nCancer){
    width=10
    p1 <- coxplot(dfm,"med",colors1)
    p2 <- coxplot(dff,"fit",colors2)
    pdf(file=outFilePath,width=width,height=4)
    suppressWarnings(grid.arrange(p1,p2, ncol=2))
    dev.off()
  }else if (nrow(dfm)>=nCancer & nrow(dff)<nCancer){
    width=5
    p1 <- coxplot(dfm,"med",colors1)
    pdf(file=outFilePath,width=width,height=4)
    suppressWarnings(print(p1))
    dev.off()
  }else if (nrow(dfm)<nCancer & nrow(dff)>=nCancer){
    width=5
    p2 <- coxplot(dff,"fit",colors2)
    pdf(file=outFilePath,width=width,height=4)
    suppressWarnings(print(p2))
    dev.off()
  }else if (nrow(dfm)<nCancer & nrow(dff)<nCancer){
    width=5
    pdf(file=outFilePath,width=width,height=4)
    temptext1 <- paste("! ",SpliceEvent,"\nNo data presented.\n","Try another one!",sep="")
    textplot(temptext1,valign="top", cex=1.3, halign= "center",col="grey")
    dev.off()
  }
}

PanNon <- function(SpliceEvent){
  outFilePath <- paste(outDir,SpliceEvent,"-",panx,".pdf",sep = "")
  width=5
  pdf(file=outFilePath,width=width,height=4)
  temptext1 <- paste("! ",SpliceEvent,"\nNo data presented.\n","Try another one!",sep="")
  textplot(temptext1,valign="top", cex=1.3, halign= "center",col="grey")
  dev.off()
}

sapply(ase,PanCox)

if (survType=="OS"){
  aseNon <- ase[!ase%in%infoAS$SpliceEvent]
  sapply(aseNon,PanNon)
}
