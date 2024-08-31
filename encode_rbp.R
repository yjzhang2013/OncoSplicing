

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

panDiff <- function(celltype, idtype, rbp){
  library(gridExtra)
  library(ggrepel)
  library(ggplot2)
  library(scales)
  
  Xtitle <- "PSI difference upon KD/KO of different RBPs"
  
  df <- readRDS(paste0(dataDir, celltype, "/", idtype, ".rds"))
  
  idmap <- unique(df$idMap)
  
  df$Type <- ifelse(df$fdr<0.05 & abs(df$psiDiff)>0.1,
                    ifelse(df$psiDiff>0.1,"Up","Down"),"Non-sig")
  df$Type <- factor(df$Type, levels = c("Up", "Non-sig", "Down"))
  df$logfdr <- ifelse(df$fdr==0, 16, df$logfdr)
  
  df <- df[order(df$logfdr*abs(df$psiDiff), decreasing = T),]
  df$text <- ifelse(seq(nrow(df))<15, df$rbp, "")
  
  
  maxdif <- max(abs(df$psiDiff))
  n = floor(-log(maxdif, 10) + 1)
  maxval <- round(maxdif,n)
  Xlimits <- c(-max(maxdif,maxval),max(maxdif,maxval))
  Xbreaks <- c(-maxval,-maxval/2,0,maxval/2,maxval)
  Xquat <- maxval/2
  
  Xtf <- (length(df$psiDiff[df$psiDiff<=0.2])>1)
  
  if (maxval>0.8 & Xtf){
    df$psiDiff <- sapply(df$psiDiff, Xtransf)
    Xlabels <- c(-maxval,-0.2,0,0.2,maxval)
    Xcolor <- c("black","red","black","red","black")
  }else{
    Xlabels <- Xbreaks
    Xcolor <- rep("black",length(Xbreaks))
  }
  
  
  ## Y-axis
  maxbh <- max(df$logfdr)
  y = floor(-log(maxbh, 10) + 1)
  Ymaxval <- round(maxbh,y)
  Ylimits <- c(0,max(maxbh,Ymaxval))
  Ybreaks <- c(0,Ymaxval*0.25,Ymaxval*0.5,Ymaxval*0.75,Ymaxval)
  Yquat <- Ymaxval*0.25
  
  Ylabels <- gsub(20, "Inf", Ybreaks)
  Ycolor <- c(rep("black",length(Ybreaks)-1), ifelse(grep("Inf", Ylabels), "red", "black"))
  cutoff <- -log(0.05,10)
  
  
  
  ggplot(df, aes(psiDiff, logfdr, fill=Type, colour = Type)) +
    ggrepel::geom_text_repel(aes(label=text), size=3, show.legend = F) +
    geom_point(shape = 21, size = 1.5, stroke=0.2, alpha=1, show.legend = T) +
    theme_bw(base_size = 10) +
    # guides(size=guide_legend(title="PSI change")) +
    # guides(shape = guide_legend(override.aes = list(size = 01))) +
    scale_fill_manual(values = c("Up"="red","Down"="blue","Non-sig"="grey70")) +
    scale_colour_manual(values = c("Up"="red","Down"="blue","Non-sig"="grey70")) +
    scale_y_continuous(limits=Ylimits,breaks=Ybreaks,labels=Ylabels,position="left") +
    scale_x_continuous(limits=Xlimits,breaks=Xbreaks,labels=Xlabels,position="bottom") +
    labs(x = Xtitle, y="-log10 ( FDR )", title = list(unique(df$idMap))) +
    # geom_text(aes(label=text), color = "grey30", check_overlap = T, size = 3)+
    geom_hline(yintercept = cutoff, size = 0.4, color="#dd7d6a",linetype = "dashed") +
    theme(title=element_text(size=8),
          axis.text.x = element_text(size=10, color=Xcolor),
          axis.text.y = element_text(size=10, color=Ycolor),
          axis.title = element_text(size=10),
          legend.title = element_text(size=10)
    )
  
  
  pdffile <- paste0(outDir, celltype, "_", idtype, ".pdf")
  ggsave(pdffile, width=6,height=5)
}




celltype <- "HepG2"
idtype <- "UBP1_ES_63866"
rbp <- "SRSF1"


dataDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/Encode/panData/"
outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/Encode/pan/"



panDiff(celltype, idtype, rbp)
