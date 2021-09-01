library(survival)
library(survminer)
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

fileDir <- "/home/webdata/spliceseq/newAdd/survival/"
outDir <- "/home/webdata/spliceseq/newAdd/ciplots/"

if (file.exists(outDir)==FALSE){
  dir.create(outDir,recursive = T)
}

infoAS <- readRDS(file="/home/webdata/spliceseq/newAdd/spliceseq_sase_med_all_survtype.rds")
infoAS$types <- paste(infoAS$cancerType,infoAS$survType,sep="_")
infoAS$infos <- paste(infoAS$cancerType,infoAS$Gene_Symbol,infoAS$SpliceEvent,infoAS$CI_Type,sep=";")

infoAS <- infoAS[order(infoAS$cancerType,infoAS$survType,infoAS$SpliceEvent),]

chunkSize <- nrow(infoAS)%/%chunkN
chunkLast <- chunkSize + nrow(infoAS)%%chunkN
infoAS$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))
infoAS <- infoAS[infoAS$chunk==chunkI,]

KMplot=function(infos){
  cat("\nchunkI-",chunkI,"-",infos)
  CancerType <- strsplit(infos,split=";")[[1]][1]
  GeneSymbol <- strsplit(infos,split=";")[[1]][2]
  SpliceEvent <- strsplit(infos,split=";")[[1]][3]
  CItype <- strsplit(infos,split=";")[[1]][4]
  outFilePath <- paste(outDir,CancerType,"-",GeneSymbol,"-",SpliceEvent,"-",CItype,".pdf",sep = "")
  
  ## survival result
  med <- round(infoX[infoX$SpliceEvent==SpliceEvent,"cutmed"],3)
  fit <- round(infoX[infoX$SpliceEvent==SpliceEvent,"cutfit"],3)
  cutoff <- list(med,fit)
  
  data <- list()
  if (!is.na(med)){
    if(!is.na(fit)){
      data[[1]] <- meddata[,c("timeY","event",SpliceEvent)]
      data[[2]] <- fitdata[,c("timeY","event",SpliceEvent)]
      ncol=2
      width=8
    }else if(is.na(fit)){
      data[[1]] <- meddata[,c("timeY","event",SpliceEvent)]
      ncol=1
      width=4
    }
    
    splots <- list()
    for (i in 1:length(data)){
      # os <- osPrev[i]
      cut <- cutoff[[i]]
      if(i==1){
        lgtitle <- "Median cutoff"
      }else{
        lgtitle <- "Optimal cutoff"
      }
      
      dataM <- data[[i]]
      
      x <- surv_fit(Surv(timeY, event) ~ dataM[,SpliceEvent], data=dataM)
      if(cut==1){
        lPSI  <- paste("PSI<",cut,"(n=", x$n[1], ")", sep = "")
        hPSI  <- paste("PSI=",cut,"(n=", x$n[2], ")", sep = "") 
      }else if(cut==0){
        lPSI  <- paste("PSI=",cut,"(n=", x$n[1], ")", sep = "")
        hPSI  <- paste("PSI>",cut,"(n=", x$n[2], ")", sep = "") 
      }else{
        lPSI  <- paste("PSI<=",cut,"(n=", x$n[1], ")", sep = "")
        hPSI  <- paste("PSI>",cut,"(n=", x$n[2], ")", sep = "") 
      }
      
      sdf <- survdiff(Surv(timeY, event) ~ dataM[,SpliceEvent], data=dataM)
      pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      pval<-format(pval, scientific = TRUE,digits = 3)
      title=paste(CancerType,"-SpliceSeq\n",SpliceEvent,"\nLog-rank p = ",pval,sep="")
      ylabName <- paste(survType,"Survival",sep=" ")
      
      splots[[i]] <- ggsurvplot(x,data=dataM,pval=F,palette=c("blue","red"),censor.shape="|", censor.size = 1.5,
                                xlab="Years",ylab=ylabName,legend.labs =c(lPSI,hPSI),
                                legend=c(0.8,1.13),legend.title = lgtitle, title=title,
                                ggtheme = theme_classic() + 
                                  theme(title=element_text(size=9),
                                        legend.title = element_text(size=10),
                                        legend.text = element_text(size=9),
                                        axis.text = element_text(size=9,color="black"),
                                        axis.title = element_text(size=9,face="bold")))
    }
    
    p <- arrange_ggsurvplots(splots,ncol=ncol,nrow=1)
    ## to avoid the production of blank first page in the output file
    ggsave(p, file = outFilePath,width=width,height=3)
    
  }else if(is.na(med)){
    library(gplots)
    width=4
    pdf(file=outFilePath,width=width,height=3)
    temptext1 <- paste("! ",SpliceEvent,"\nNo data presented.\n","Try another one!",sep="")
    textplot(temptext1,valign="top", cex=1.3, halign= "center",col="grey")
    dev.off()
  }
}

types <- infoAS$types[!duplicated(infoAS$types)]

for (tt in types){
  cancerType <- strsplit(tt,split="_")[[1]][1]
  survType <- strsplit(tt,split="_")[[1]][2]
  cat("\nPloting:",cancerType,survType)
  
  infoX <- infoAS[infoAS$cancerType==cancerType & infoAS$survType==survType,]
  load(paste(fileDir,"spliceseq_med_data_",survType,"_",cancerType,".Rdata",sep = ""))
  meddata <- mydata
  load(paste(fileDir,"spliceseq_fit_data_",survType,"_",cancerType,".Rdata",sep = ""))
  fitdata <- mydata
  
  infos <- infoX$infos
  sapply(infos,KMplot)
}
