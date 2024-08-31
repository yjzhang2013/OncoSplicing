encode_plot <- function(celltype, idtype, rbp){
  res <- readRDS(paste0(encDir, "plotData/", idtype, ".rds"))
  cellRbp <- paste(celltype, rbp, sep="_")
  xTitle <- paste0("Splicing upon target RBP ", rbp, " in cell ", celltype)
  output <- paste0(encDir, "plot/", celltype,"_", idtype, "_", rbp, ".pdf")
  source(paste0(encDir, "script/encode_vis.R"))
  source(paste0(encDir, "script/geom_arch.R"))
  
  
  gtf_df <- res[["gtf_df"]]
  region_df <- res[["region_df"]][[cellRbp]]
  bed_df <- res[["bed_df"]][[cellRbp]]
  
  filetype <- gsub(".*_(.*)_.*", "\\1", region_df$fileName)
  filename <- gsub("_(.*)_", "_", region_df$fileName)
  region_df$fileName <- paste0(filetype, "_", filename)
  region_df <- region_df[order(region_df$fileName),]
  
  sj <- region_df[region_df$track_type=="junction",]
  sjs <- length(unique(paste(sj$seqnames, sj$start, sj$end, sep=":")))
  
  if (nrow(sj)!=4*sjs){
    sjcheck <- sj %>% filter(!duplicated(paste(seqnames, start, end, sep=":")))
    
    sjcheck <- rbind(sjcheck, sjcheck, sjcheck, sjcheck)
    sjcheck <- sjcheck %>% mutate(score=0, fileName=rep(unique(sj$fileName), sjs))
    sjcheck <- sjcheck %>% filter(!paste0(seqnames, start, end, fileName)%in%paste0(sj$seqnames, sj$start, sj$end, sj$fileName))
    region_df <- rbind(region_df, sjcheck)
  }
  
  
  Input_gene <- unique(gtf_df$gene_name)
  
  if(!is.null(bed_df)){
    w <- 11
    h <- 5
  }else{
    w <- 8
    h <- 4
  }
  
  options(scipen = 7)
  s <- Sys.time()
  trackVisProMax2(Input_gtf = gtf_df,
                  bed.df = bed_df,
                  region.df = region_df,
                  Input_gene = Input_gene,
                  xTitle = xTitle, 
  )
  ggplot2::ggsave(output, width = w, height = h)
  e <- Sys.time()
  e-s
}

encDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/Encode/"


idtype <- "RBM39_ES_59250"
celltype <- "HepG2"
rbp <- "RBM39"

encode_plot(celltype, idtype, rbp)


trans_fill_col = RColorBrewer::brewer.pal(n =3, name = "Paired")

