mapas_plot <- function(idtype, rbp){
  res <- readRDS(paste0(mapDir, "plotData2/", idtype, ".rds"))
  xTitle <- paste0("Motifs/peaks of RBP ", rbp, "")
  output <- paste0(mapDir, "plot/", idtype, "_", rbp, ".pdf")
  source(paste0(mapDir, "script/encode_vis.R"))
  source(paste0(mapDir, "script/geom_arch.R"))
  
  region_df <- NULL
  gtf_df <- res[["gtf_df"]]
  bed_df <- res[["bed_df"]][[rbp]]
  
  Input_gene <- unique(gtf_df$gene_name)
  
  
  
  options(scipen = 7)
  s <- Sys.time()
  trackVisProMax2(Input_gtf = gtf_df,
                  bed.df = bed_df,
                  region.df = region_df,
                  Input_gene = Input_gene,
                  xTitle = xTitle, 
  )
  ggplot2::ggsave(output, width = 10, height = 4)
  e <- Sys.time()
  e-s
}


mapDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/MapAS/"


idtype <- "alt_5prime_178640"
rbp <- "ILF2"

mapas_plot(idtype, rbp)

