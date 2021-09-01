library(survival)
library(survminer)
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
outDir <- paste("/home/webdata/spliceseq/newplots/",CancerType,"/",sep="")

if (file.exists(outDir)==FALSE){
  dir.create(outDir,recursive = T)
}

infoAS <- readRDS(file=paste(fileDir,"spliceseq_info_cutoff_",CancerType,".rds",sep=""))
chunkSize <- nrow(infoAS)%/%chunkN
chunkLast <- chunkSize + nrow(infoAS)%%chunkN
infoAS$chunk <- c(rep(1:(chunkN-1),each=chunkSize),rep(chunkN,chunkLast))
infoAS <- infoAS[infoAS$chunk==chunkI,]


cTypePFS <- c("PCPG","PRAD","TGCT","THCA","THYM")

if (CancerType %in% cTypePFS){
  ylabName <- "PFI Survival"
}else{
  ylabName <- "OS Survival"
}


KMplot=function(SpliceEvent){
  outFilePath <- paste(outDir,CancerType,"-",SpliceEvent,"-KMplot.pdf",sep = "")
  
  ## survival result
  med <- round(infoAS[infoAS$SpliceEvent==SpliceEvent,"cutmed"],3)
  fit <- round(infoAS[infoAS$SpliceEvent==SpliceEvent,"cutfit"],3)
  cutoff <- list(med,fit)
  
  if (!is.na(med)){
    if(!is.na(fit)){
      osPrev <- c("spliceseq_med_data_os_","spliceseq_fit_data_os_")
      ncol=2
      width=8
    }else if(is.na(fit)){
      osPrev <- c("spliceseq_med_data_os_")
      ncol=1
      width=4
    }
    
    splots <- list()
    for (i in 1:length(osPrev)){
      os <- osPrev[i]
      cut <- cutoff[[i]]
      if(i==1){
        lgtitle <- "Median cutoff"
      }else{
        lgtitle <- "Optimal cutoff"
      }
      
      load(paste(fileDir,os,CancerType,".Rdata",sep = ""))
      dataM <- mydata[,c("timeY","event",SpliceEvent)]
      
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
      
      splots[[os]] <- ggsurvplot(x,data=dataM,pval=F,palette=c("blue","red"),censor.shape="|", censor.size = 1.5,
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
  
  return(outFilePath)
}

ase <- infoAS$SpliceEvent
sapply(ase,KMplot)
