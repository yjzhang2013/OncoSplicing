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

fileDir <- "/home/webdata/spladder/dataplot/"
outDir <- paste("/home/webdata/spladder/tnplots/",CancerType,"/",sep="")

if (file.exists(outDir)==FALSE){
  dir.create(outDir,recursive = T)
}

psi <- readRDS(paste(fileDir,"spladder_data_psi_",CancerType,".rds",sep=""))
infoAS <- readRDS(file=paste(fileDir,"spladder_info_cutoff_",CancerType,".rds",sep=""))

chunkSize <- nrow(infoAS)%/%chunkN
chunkLast <- chunkSize + nrow(infoAS)%%chunkN

infoAS$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))
infoAS <- infoAS[infoAS$chunk==chunkI,]


TNplot=function(SpliceEvent){
  cat("\nchunkI",chunkI,"-",SpliceEvent)
  idType <- infoAS[infoAS$SpliceEvent==SpliceEvent,"idType"]
  outFilePath <- paste(outDir,CancerType,"-",SpliceEvent,"-TNplot.pdf",sep = "")  

  df <- psi[,c("tissueType",SpliceEvent)]
  df <- na.omit(df)
  tb <- data.frame(table(df$tissueType))
  names(df)[2] <- "PSI"
  
  t <- paste("TCGA",CancerType,"T",sep="-")
  t <- as.character(tb$Var1[grepl(t,tb$Var1)])
  n <- paste("TCGA",CancerType,"N",sep="-")
  n <- as.character(tb$Var1[grepl(n,tb$Var1)])
  g <- as.character(tb$Var1[grepl("GTEx",tb$Var1)])
  
  if(nrow(tb)==3){
    width=4
    dfn <- df[df$tissueType==n,"PSI"]
    dft <- df[df$tissueType==t,"PSI"]
    dfg <- df[df$tissueType==g,"PSI"]
    x <- tryCatch(wilcox.test(dft,dfn)$p.value,error = function(e) return(NA))
    y <- tryCatch(wilcox.test(dft,dfg)$p.value,error = function(e) return(NA))
    x <- format(x, scientific = TRUE,digits = 3)
    y <- format(y, scientific = TRUE,digits = 3)
    title=paste(CancerType,"-SplAdder\n",idType,"\nTumor-Normal, p = ",x,"\nTumor-GTEx, p = ",y,sep="")
    gtex <- paste(g,"\nn=",length(dfg),sep="")
    norm <- paste(n,"\nn=",length(dfn),sep="")
    tumor <- paste(t,"\nn=",length(dft),sep="")
    Xlimits <- c(g,n,t)
    Xlabels <- c(gtex,norm,tumor)
    color <- c("grey80","lightblue","#FF6A6A")
  }else if(nrow(tb)==2){
    width=3
    if(length(n)==0){
      dft <- df[df$tissueType==t,"PSI"]
      dfg <- df[df$tissueType==g,"PSI"]
      y <- tryCatch(wilcox.test(dft,dfg)$p.value,error = function(e) return(NA))
      y <- format(y, scientific = TRUE,digits = 3)
      title=paste(CancerType,"-SplAdder\n",idType,"\nTumor-GTEx, p = ",y,sep="")
      gtex <- paste(g,"\nn=",length(dfg),sep="")
      tumor <- paste(t,"\nn=",length(dft),sep="")
      Xlimits <- c(g,t)
      Xlabels <- c(gtex,tumor)
      color <- c("grey80","#FF6A6A")
    }else if(length(g)==0){
      dfn <- df[df$tissueType==n,"PSI"]
      dft <- df[df$tissueType==t,"PSI"]
      x <- tryCatch(wilcox.test(dft,dfn)$p.value,error = function(e) return(NA))
      x <- format(x, scientific = TRUE,digits = 3)
      title=paste(CancerType,"-SplAdder\n",idType,"\nTumor-Normal, p = ",x,sep="")
      norm <- paste(n,"\nn=",length(dfn),sep="")
      tumor <- paste(t,"\nn=",length(dft),sep="")
      Xlimits <- c(n,t)
      Xlabels <- c(norm,tumor)
      color <- c("lightblue","#FF6A6A")
    }else if(length(t)==0){
      dfn <- df[df$tissueType==n,"PSI"]
      dfg <- df[df$tissueType==g,"PSI"]
      title=paste(CancerType,"-SplAdder\n",idType,sep="")
      norm <- paste(n,"\nn=",length(dfn),sep="")
      gtex <- paste(g,"\nn=",length(dfg),sep="")
      Xlimits <- c(n,g)
      Xlabels <- c(norm,gtex)
      color <- c("lightblue","grey80")
    }
  }else if(nrow(tb)==1){
    width=3
    title=paste(CancerType,"-SplAdder\n",idType,sep="")
    if(length(t)==1){
      dft <- df[df$tissueType==t,"PSI"]
      tumor <- paste(t,"\nn=",length(dft),sep="")
      Xlimits <- c(t)
      Xlabels <- c(tumor)
      color <- c("#FF6A6A")
    }else if(length(g)==1){
      dfg <- df[df$tissueType==g,"PSI"]
      gtex <- paste(g,"\nn=",length(dfg),sep="")
      Xlimits <- c(g)
      Xlabels <- c(gtex)
      color <- c("grey80")
    }else if(length(n)==1){
      dfn <- df[df$tissueType==n,"PSI"]
      norm <- paste(n,"\nn=",length(dfn),sep="")
      Xlimits <- c(n)
      Xlabels <- c(norm)
      color <- c("lightblue")
    }
  }
  
  pdf(file=outFilePath,width=width,height=3.7)
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
  suppressWarnings(print(p))
  suppressMessages(dev.off())
}

ase <- infoAS$SpliceEvent
sapply(ase,TNplot)
