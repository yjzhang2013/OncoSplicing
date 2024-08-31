corrTest <- function(cancerType, ase, sf){
  count <- sf[[cancerType]]
  psi <- ase[[cancerType]]
  
  cor <- psych::corr.test(psi, count, use="pairwise", adjust="fdr")
  
  # options(digits = 4, scipen = 0)
  df <- data.frame(Cancer_Type = cancerType,
                   All_Samples=length(count),
                   Samples = cor$n,
                   PCT = round(cor$n/length(count), 2),
                   Correlation = round(cor$r, 3),
                   P_Value = formattable::scientific(cor$p, digits = 2),
                   FDR = formattable::scientific(cor$p.adj, digits = 2)
  )
  
  return(df)
}


corrTable <- function(omicsA,omicsB, inputA, inputB){
  ase <- readRDS(paste0(fileDir, tolower(omicsA), "/", inputA, ".rds"))
  sf <- readRDS(paste0(fileDir, tolower(omicsB), "/", inputB, ".rds"))
  cancerTypes <- names(sf)[names(sf)%in%names(ase)]
  res <- lapply(cancerTypes, corrTest, sf, ase)
  res <- data.frame(do.call(rbind, res))
  tableName <- paste0(outDir, omicsA,"_", omicsB, "_", inputA, "_", inputB, ".csv")
  
  write.csv(res, tableName, row.names = F)
}


#
fileDir <- "/mnt/webData/CoExp/"
outDir <- "/mnt/webData/CoExp/csv/"

fileDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/CoExp/"
outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/CoExp/csv/"

# omicsA <- "mRNA"
# omicsB <- "SpliceSeq"
# inputA <- "A1CF"
# inputB <- "A2M_ES_20222"
# corrTable(omicsA, omicsB, inputA,  inputB)

corrTable("mRNA", "SpliceSeq", "A1CF",  "A2M_ES_20222")
