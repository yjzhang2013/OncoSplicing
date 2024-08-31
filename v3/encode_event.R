

## def function for x-axis transform
Xtransf <- function(data, maxval){
  xCut <- 0.2
  Xquat <- maxval/2
  
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
  library(dplyr)
  
  Xtitle <- paste0("PSI difference upon KD/KO of ", rbp, " in ", celltype)
  idmap <- rbp
  
  df <- readRDS(paste0(dataDir, celltype, "/", rbp, ".rds"))
  df$Type <- ifelse(df$fdr<0.05 & abs(df$psiDiff)>0.1,
                    ifelse(df$psiDiff>0.1,"Up","Down"),"Non-sig")
  
  colors <- c("red", "grey70", "blue")[c("Up"%in%unique(df$Type), "Non-sig"%in%unique(df$Type), "Down"%in%unique(df$Type))]
  
  df <- df %>% 
    group_by(Type) %>% 
    mutate(Type=ifelse(Type!="Non-sig", paste0(Type, " (",n(),")"), Type))
  
  df$Type <- factor(df$Type, levels = rev(unique(df$Type)[order(unique(df$Type))]))
  
  
  df$logfdr <- ifelse(df$fdr==0, 18, df$logfdr)
  
  df <- df[order(df$logfdr*abs(df$psiDiff), decreasing = T),]
  df$text <- ifelse(seq(nrow(df))<15, df$idMap, "")
  
  
  maxdif <- max(abs(df$psiDiff))
  n = floor(-log(maxdif, 10) + 1)
  maxval <- round(maxdif, n)
  Xlimits <- c(-max(maxdif,maxval),max(maxdif,maxval))
  Xbreaks <- c(-maxval,-maxval/2,0,maxval/2,maxval)
  
  
  Xtf <- (length(df$psiDiff[df$psiDiff<=0.2])>1)
  
  if (maxval>0.8 & Xtf & F){
    df$psiDiff <- sapply(df$psiDiff, Xtransf, maxval)
    Xlabels <- c(-maxval, -0.2, 0, 0.2, maxval)
    Xcolor <- c("black", "red", "black", "red", "black")
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
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = colors) +
    scale_y_continuous(limits=Ylimits,breaks=Ybreaks,labels=Ylabels,position="left") +
    scale_x_continuous(limits=Xlimits,breaks=Xbreaks,labels=Xlabels,position="bottom") +
    labs(x = Xtitle, y="-log10 ( FDR )", title = idmap) +
    # geom_text(aes(label=text), color = "grey30", check_overlap = T, size = 3)+
    geom_hline(yintercept = cutoff, size = 0.4, color="#dd7d6a",linetype = "dashed") +
    theme(title=element_text(size=8),
          axis.text.x = element_text(size=10, color=Xcolor),
          axis.text.y = element_text(size=10, color=Ycolor),
          axis.title = element_text(size=10),
          legend.title = element_text(size=10)
    )
  
  
  pdffile <- paste0(outDir, celltype, "_", rbp, ".pdf")
  ggsave(pdffile, width=6,height=5)
}



celltype <- "K562"
idtype <- "exon_skip_497057"
rbp <- "SRSF1"


dataDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/Encode/panData/"
outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/Encode/pan/"



panDiff(celltype, idtype, rbp)

