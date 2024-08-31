plot <- function(omicsA, omicsB, cancerType, inputA, inputB){
  ase <- readRDS(paste0(fileDir, tolower(omicsA), "/", inputA, ".rds"))
  sf <- readRDS(paste0(fileDir, tolower(omicsB), "/", inputB, ".rds"))
  
  count <- sf[[cancerType]]
  psi <- ase[[cancerType]]
  
  xlab <- ifelse(omicsB=="mRNA", paste0(inputB, " (log2(TPM+1))"), paste0(inputB, " (PSI)"))
  ylab <- ifelse(omicsA=="mRNA", paste0(inputA, " (log2(TPM+1))"), paste0(inputA, " (PSI)"))
  title <- paste0(cancerType, "(", omicsA, " vs ", omicsB, ")")
  
  options(scipen = -3)
  p <-
    ggpubr::ggscatter(data.frame(count, psi), x = "count", y = "psi",
                      add = "reg.line", conf.int = TRUE,
                      add.params = list(fill = "Salmon", color="red", size=0.6),
                      color = "DodgerBlue",
                      size = 0.8,
                      title = title,
                      xlab = xlab,
                      ylab = ylab,
                      family = "Arial") +
    # scale_y_continuous(limits = c(0,1.1), breaks = c(0,1,0.5)) +
    scale_x_continuous(labels = scales::label_number()) +
    scale_y_continuous(labels = scales::label_number()) +
    ggpubr::stat_cor(method = "pearson", color="red", digits = 4,
                     label.x = min(count, na.rm = T), label.y = min(psi, na.rm = T))
  
  pdfName <- paste0(outDir, omicsA, "_", omicsB, "_", cancerType, "_", inputA, "_", inputB, ".pdf")
  #suppressWarnings(suppressMessages(ggplot2::ggsave(pdfName, width = 4, height = 4)))
  ggplot2::ggsave(pdfName, width = 4, height = 4)
}


# fileDir <- "/mnt/webData/CoExp/"
# outDir <- "/mnt/webData/CoExp/pdf/"

fileDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/CoExp/"
outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/CoExp/pdf/"

omicsA <- "mRNA"
omicsB <- "SpliceSeq"
inputA <- "ESRP2"
inputB <- "PLEKHM2_ES_767"
cancerType <- "TCGA-KIRC"
plot(omicsA, omicsB, cancerType,inputA, inputB)


