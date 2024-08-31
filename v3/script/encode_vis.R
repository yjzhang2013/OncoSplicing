library(ggplot2)
suppressMessages(library(dplyr))


globalVariables(c("Freq","dist", "element_line", "exon_len", "facetted_pos_scales", "geom_label",
                  "ggplotGrob", "guide_legend", "label_pos", "labeller",
                  "margin", "sampleName", "scale_color_gradientn",
                  "scale_color_manual", "scale_fill_gradientn", "scale_x_continuous","smin_new",
                  "scale_y_continuous", "score", "smax", "smax_label", "smax_new", "smin",
                  "theme_void", "track_type", "xend", "yend", "ymax", "ymin"))


trackVisProMax2 <- function(Input_gtf = NULL,
                            region.df = NULL,
                            bed.df = NULL,
                            Input_gene = NULL,
                            Loop_curve_geom = "geom_arch",
                            xTitle = "",
                            signal_layer_bw_params = list(),
                            peaks_layer_params = list(),
                            junc_layer_combined = TRUE,
                            signal_layer_junction_params = list(),
                            signal_layer_junction_label_params = list(),
                            reverse_y_vars = NULL,
                            collapse_trans = FALSE,
                            exon_width = 0.8,
                            peak_width = 0.5,
                            arrow_rel_len_params_list = list(),
                            trans_exon_arrow_params = list(),
                            trans_exon_col_params = list(),
                            gene_dist_mark_params = list(),
                            gene_dist_mark_text_params = list(),
                            signal_range_pos = c(0.1, 0.85),
                            panel_size_setting = list(),
                            fixed_column_range = TRUE,
                            signal_range_label_params = list(),
                            by_layer_x = FALSE,
                            by_layer_y = FALSE,
                            column_strip_setting_list = list(),
                            row_strip_setting_list = list(),
                            column_strip_text_setting_list = list(),
                            row_strip_text_setting_list = list(),
                            add_gene_region_label = FALSE,
                            base_size = 14,
                            panel.spacing = c(0.2,0),
                            sample_fill_col = c("#D51F26","#D51F26", "#272E6A", "#272E6A"),
                            # peak_fill_col = c("#208A42","#89288F","#F47D2B"),
                            trans_fill_col = NULL,
                            Intron_line_type = "line"
){
  options(warn=-1)
  options(dplyr.summarise.inform = FALSE)
  
  gtf <- Input_gtf
  gene = Input_gene
  
  # ==============================================================================
  # 2_add trans and chromo facet
  # ==============================================================================
  add_facet_name = c("Peaks/motif","Splice Event")
  tran_facet <- data.frame(seqnames = NA,start = NA,end = NA,score = NA,
                           fileName = rep(add_facet_name, each = 1),
                           track_type = NA,
                           gene = gene)
  
  
  if (!is.null(region.df)){
    # merge
    tmp2 <- plyr::rbind.fill(region.df, tran_facet) %>% unique()
    
    # gene orders
    tmp2$gene <- factor(tmp2$gene, levels = unique(region.df$gene))
    
    # sample orders
    tmp2$fileName <- factor(tmp2$fileName,
                            levels = unique(c(unique(region.df$fileName), add_facet_name)))
  }else{
    tmp2 <- tran_facet
  }
  
  # x = 1
  if(!is.null(bed.df)){
    bed.df.new <- bed.df %>%
      mutate(fileName = "Peaks/motif",
             gene = gene)
  }
  
  
  # ==============================================================================
  # 4_gene structures
  # ==============================================================================
  tmp_gtf <- gtf %>%
    dplyr::filter(gene_name %in% Input_gene & type != "gene") %>%
    mutate(gene = gene_name)
  
  transcript.df <- tmp_gtf %>%
    mutate(y = 1,
           ymin = 1 - exon_width*0.5,
           ymax = 1 + exon_width*0.5,
           gene = gene_name, fileName = "Splice Event")
  
  
  # =============================================
  # transcript arrow
  # =============================================
  # x = 1
  tmp <- transcript.df %>%
    dplyr::filter(gene_id == Input_gene & type == "transcript")
  
  # create segment
  plyr::ldply(1:nrow(tmp),function(x){
    tmp1 <- tmp[x,]
    ypos = unique(transcript.df[which(transcript.df$transcript_id == tmp1$transcript_id),
                                c("y")])
    
    if(Intron_line_type == "chevron"){
      tmp_exon <- gtf %>%
        filter(transcript_id ==  tmp1$transcript_id & type == "exon") %>%
        arrange(start,end)
      
      xstart = tmp_exon$end[1:nrow(tmp_exon) - 1]
      xend = tmp_exon$start[2:nrow(tmp_exon)]
      seg_data <- data.frame(x = c(xstart,(xstart + xend)/2),
                             xend = c((xstart + xend)/2,xend),
                             y = rep(c(ypos,ypos + exon_width*0.25),each = nrow(tmp_exon) - 1),
                             yend = rev(rep(c(ypos,ypos + exon_width*0.25),each = nrow(tmp_exon) - 1)),
                             transcript_id = tmp1$transcript_id
      )
    }else if(Intron_line_type == "line"){
      seg_data <-
        do.call(createSegment,modifyList(list(xPos = c(tmp1$start,tmp1$end),
                                              yPos = rep(ypos,2),
                                              rel_len = 0.08),
                                         arrow_rel_len_params_list)) %>%
        mutate(transcript_id = tmp1$transcript_id)
    }
    
    # get group info
    seg_data$gene <- tmp1$gene
    seg_data$ends <- ifelse(unique(tmp$strand) == "+","last","first")
    seg_data$fileName <- "Splice Event"
    
    return(seg_data)
  }) -> final_arrow_data
  
  
  
  # trans_exon_arrow_params = list(length = 1,fill = "grey60",color = "grey60",linewidth = 0.75)
  trans_arrow_layer <- lapply(unique(final_arrow_data$transcript_id), function(x){
    tmp <- final_arrow_data[which(final_arrow_data$transcript_id == x),]
    
    # arrow layer
    if(Intron_line_type == "chevron"){
      do.call(geom_segment,
              modifyList(list(data = tmp,
                              aes(x = x,xend = xend,
                                  y = y,yend = yend),
                              linewidth = 0.75,
                              color = "grey60"),trans_exon_arrow_params))
    }else{
      do.call(geom_segment,
              modifyList(list(data = tmp,
                              aes(x = x,xend = xend,
                                  y = y,yend = yend),
                              linewidth = 0.75,
                              arrow = arrow(type = "closed",
                                            length = unit(1,"mm"),
                                            ends = unique(tmp$ends)),
                              arrow.fill = "grey60",
                              color = "grey60"),trans_exon_arrow_params))
    }
  })
  
  
  if(!is.null(trans_exon_col_params$mapping)){
    trans_mapping <- list(data = transcript.df %>% dplyr::filter(type != "transcript"),
                          mapping= aes(xmin = start,xmax = end,
                                       ymin = ymin, ymax = ymax,
                                       fill = "orange"),
                          color = "grey60")
  }else{
    trans_mapping <- list(data = transcript.df %>% dplyr::filter(type != "transcript"),
                          mapping= aes(xmin = start,xmax = end,
                                       ymin = ymin, ymax = ymax),
                          fill = "orange",
                          color = "grey60")
  }
  
  
  trans_struct_layer <- do.call(geom_rect,
                                modifyList(trans_mapping,trans_exon_col_params))
  
  # ==============================================================================
  # 5_segment and arrow data for chromosome label and region length
  # ==============================================================================
  segment.df <- transcript.df %>%
    group_by(fileName,gene,seqnames,strand)
  
  # add segment info
  segment.df <- segment.df %>%
    summarise(start = min(start),end = max(end))
  
  # two segment lines position
  segment.df <- segment.df %>%
    mutate(ar1_end = start + (end - start)/3,
           ar2_start = end - (end - start)/3)
  
  # arrow data on gene structures
  # x = 1
  plyr::ldply(1:nrow(segment.df),function(x){
    tmp <- segment.df[x,] %>%
      # add chr prefix
      mutate(seqnames = if_else(startsWith(as.character(seqnames),"chr"),
                                seqnames,paste("chr",seqnames,sep = "")))
    # calculate arrow positions
    # t_num = table(data.frame(unique(transcript.df[,c("gene","transcript_id")]))$gene)
    t_num = table(data.frame(unique(final_arrow_data[,c("gene","y")]))$gene)
    
    res <- data.frame(start = c(tmp$start,tmp$end - (tmp$end - tmp$start)/3),
                      end = c(tmp$start + (tmp$end - tmp$start)/3,tmp$end)) %>%
      mutate(fileName = tmp$fileName,gene = tmp$gene,strand = tmp$strand,
             .before = start) %>%
      mutate(label = paste(tmp$seqnames,": ",round((tmp$end - tmp$start)/10^3,digits = 2),
                           " kb",sep = ""),
             label_pos = (tmp$end + tmp$start)/2) %>%
      mutate(y = t_num[tmp$gene] + 1)
    
    # add arrow direction
    res <- res %>%
      mutate(ends = if_else(strand == "+","last","first"))
    
    return(res)
  }) -> arrow.df
  
  seg_arrow_layer <- lapply(unique(arrow.df$gene), function(x){
    tmp <- arrow.df %>% dplyr::filter(gene == x)
    seg_ar <-
      do.call(geom_segment,
              modifyList(list(data = tmp,
                              aes(x = start,xend = end,
                                  y = if(collapse_trans == TRUE){y = 2}else{y},
                                  yend = if(collapse_trans == TRUE){y = 2}else{y}
                              ),
                              linewidth = 0.3,
                              color = "black",
                              arrow = arrow(ends = tmp$ends,
                                            length = unit(2,"mm"))),
                         gene_dist_mark_params))
    return(seg_ar)
  })
  
  text_layer <- do.call(geom_text,
                        modifyList(list(data = arrow.df,
                                        aes(x = label_pos,
                                            y = if(collapse_trans == TRUE){y = 2}else{y},
                                            label = label),
                                        color = "black",
                                        size = 3),
                                   gene_dist_mark_text_params))
  
  # ==============================================================================
  # 6_text and signal ranges for each panel
  # ==============================================================================
  if(!is.null(region.df)){
    # signal_range_pos = c(0.9,0.9)
    rg_xpos <- tmp2 %>%
      dplyr::filter(!(fileName %in% add_facet_name)) %>%
      group_by(gene) %>%
      summarise(xmin = min(start),xmax = max(end)) %>%
      ungroup() %>% mutate(rg_xpos = (xmax - xmin)*signal_range_pos[1] + xmin) %>%
      dplyr::select(gene,rg_xpos)
    
    # merge with rg_xpos
    rg_info <- tmp2 %>%
      dplyr::filter(!(fileName %in% add_facet_name)) %>%
      group_by(fileName,gene,track_type) %>%
      summarise(smin = min(score),smax = max(score)) %>%
      ungroup() %>%
      # mutate(rg_ypos = (smax - smin)*signal_range_pos[2] + smin) %>%
      left_join(.,rg_xpos,by = "gene")
    
    # whether dplyr::filter junction
    if(junc_layer_combined == TRUE){
      rg_info <- rg_info[which(rg_info$track_type != "junction"),]
    }
    
    
    # add new range column
    # x = 1
    if(fixed_column_range == TRUE){
      plyr::ldply(1:length(unique(rg_info$gene)),function(x){
        tmp <- rg_info %>% dplyr::filter(gene == unique(rg_info$gene)[x]) %>%
          mutate(smax = max(smax),smin = min(smin))
      }) -> new_range
    }else{
      new_range <- rg_info
    }
    
    # add label yPos
    new_range <- new_range %>% mutate(rg_ypos = (smax - smin)*signal_range_pos[2] + smin)
    
    
    # add range text label column
    new_range <- new_range %>%
      mutate(smax_value = ceiling(smax),
             smin_value = floor(smin),
             smax_label = paste("[",as.integer(floor(smin)),"-",as.integer(ceiling(smax)),"]",sep = ""))
    
    # ============================
    # whether reverse y axis
    new_range <- new_range %>%
      mutate(yr_type = ifelse(fileName %in% reverse_y_vars,"reverse","identity"))
    
    # order
    new_range$fileName <- factor(new_range$fileName,levels = levels(tmp2$fileName))
    
    # layer
    # signal_range_label_params = list(color = "black",size = 4)
    range_label_layer <-
      do.call(zplyr::geom_abs_text,
              modifyList(list(data = new_range,
                              aes(xpos = signal_range_pos[1],
                                  ypos = signal_range_pos[2],
                                  label = smax_label),
                              color = "black",
                              size = 4),signal_range_label_params))
    
    # mutate panel order
    new_range <- new_range %>% arrange(fileName,
                                       gene) %>%
      mutate(panel_num = 1:nrow(new_range))
    
    # trans panel y range
    tarns_pos_y_range <- unique(arrow.df[,c("gene","y")])
    
    # panel_range_layer
    iter_loop <- 1:(nrow(new_range) + length(add_facet_name)*(length(unique(tmp2$gene))))
    
    panel_range_layer <- lapply(iter_loop, function(x){
      if(x <= nrow(new_range)){
        tmp <- new_range %>% dplyr::filter(panel_num == x)
        if(tmp$yr_type == "reverse"){
          sy <- scale_y_continuous(limits = c(tmp$smax_value,tmp$smin_value),trans = "reverse")
        }else{
          sy <- scale_y_continuous(limits = c(tmp$smin_value,tmp$smax_value))
        }
      }else{
        if(!is.null(bed.df) & x %in%
           c((nrow(new_range) + 1):(nrow(new_range) + length(unique(bed.df$sampleName))))){
          sy <- scale_y_continuous(limits = c(0,length(unique(bed.df$sampleName)) + 1))
        }else{
          if(collapse_trans == TRUE){
            if(!is.null(Input_gene)){
              sy <- scale_y_continuous(limits = c(0,1 + 1.5))
            }else{
              sy <- scale_y_continuous(limits = c(0,2.5))
            }
          }else{
            if(!is.null(Input_gene)){
              if(!is.null(Input_bed)){
                index <- x - nrow(new_range) - length(unique(bed.df$sampleName))
              }else{
                index <- x - nrow(new_range)
              }
              sy <- scale_y_continuous(limits = c(0,tarns_pos_y_range[index,]$y + 0.5))
            }else{
              tmp_p <- final_arrow_data %>% dplyr::select(gene,y) %>%
                group_by(gene) %>% summarise(y = max(y))
              sy <- scale_y_continuous(limits = c(0,tmp_p$y[x - nrow(new_range)] + 1.5))
              
            }
          }
        }
      }
      return(sy)
    })
  }else{
    # ==========================================================================
    range_label_layer <- NULL
    
    # trans panel y range
    tarns_pos_y_range <- unique(arrow.df[,c("gene","y")])
    
    iter_loop <- 1:(length(add_facet_name)*(length(unique(transcript.df$gene))))
    panel_range_layer <- lapply(iter_loop, function(x){
      if(collapse_trans == TRUE){
        sy <- scale_y_continuous(limits = c(0,2.5))
      }else{
        sy <- scale_y_continuous(limits = c(0,tarns_pos_y_range[x,]$y + 1.5))
      }
      return(sy)
    })
  }
  
  
  panel_x_range_layer <- NULL
  
  # ==============================================================================
  # 10_facet settings
  # ==============================================================================
  facet_strips <- ggh4x::strip_themed(
    background_x = do.call(ggh4x::elem_list_rect,
                           modifyList(list(colour = rep("white",2*length(unique(region.df$gene)))),
                                      column_strip_setting_list)),
    by_layer_x = by_layer_x,
    text_x = do.call(ggh4x::elem_list_text,
                     modifyList(list(),
                                column_strip_text_setting_list)),
    
    
    background_y = do.call(ggh4x::elem_list_rect,
                           modifyList(list(colour = rep("white",2*length(unique(region.df$gene)))),
                                      row_strip_setting_list)),
    by_layer_y = by_layer_y,
    text_y = do.call(ggh4x::elem_list_text,
                     modifyList(list(),
                                row_strip_text_setting_list)),
  )
  
  # whether supply with parames
  if(length(column_strip_setting_list) == 0 &
     length(column_strip_text_setting_list) == 0 &
     length(row_strip_setting_list) == 0 &
     length(row_strip_text_setting_list) == 0){
    facet_strips <- ggh4x::strip_nested()
  }else{
    facet_strips <- facet_strips
  }
  
  facet_var <- fileName~gene
  
  new_column_label = NULL
  
  facet_layer <- do.call(ggh4x::facet_nested,
                         modifyList(list(facet_var,
                                         scales = "free",
                                         nest_line = element_line(linetype = "solid"),
                                         independent = "y",
                                         switch = "y",
                                         strip = facet_strips,
                                         labeller = labeller(gene = new_column_label)),
                                    list()))
  
  # width and height of each panel
  # panel_size_setting = list()
  if(length(panel_size_setting) == 0){
    panel_size_layer <- NULL
  }else{
    panel_size_layer <- do.call(ggh4x::force_panelsizes,modifyList(
      list(respect = TRUE),panel_size_setting))
  }
  
  # ==============================================================================
  # 11_peaks facet plot
  # ==============================================================================
  # peaks_layer_params = list()
  if(!is.null(bed.df)){
    peaks_layer <- do.call(geom_rect,
                           modifyList(list(data = bed.df.new,
                                           aes(xmin = start,xmax = end,
                                               ymin = ymin,ymax = ymax,
                                               fill = sampleName),
                                           color = NA),
                                      peaks_layer_params))
  }else{
    peaks_layer <- NULL
  }
  
  # ==============================================================================
  # 12_signal layer
  # ==============================================================================
  # bigwig layer
  # bw_geom = "rect"
  if(!is.null(region.df)){
    bw_data <- tmp2[which(tmp2$track_type == "bigwig"),]
    bw_data$start <- c(bw_data$start[1],bw_data$start[2:nrow(bw_data)]-1)
    
    signal_layer_geom_rect <- do.call(geom_rect,
                                      modifyList(
                                        list(data = bw_data,
                                             mapping = aes(xmin = start,xmax = end,
                                                           ymin = 0, ymax = score,
                                                           fill = fileName),
                                             color = NA,size = 0,
                                             show.legend = FALSE),
                                        signal_layer_bw_params))
  }else{
    signal_layer_geom_rect <- NULL
  }
  
  
  # ============================================================================
  # junction data process
  if(!is.null(region.df)){
    # junc_layer_combined = FALSE
    range_bw <- new_range[which(new_range$track_type == "bigwig"),]
    junc_data <- tmp2[which(tmp2$track_type == "junction"),] %>%
      mutate(dist = end - start)
    
    # g = 1;f = 1
    file_name <- unique(junc_data$fileName)
    gene <- unique(junc_data$gene)
    plyr::ldply(seq_along(gene),function(g){
      tmp_g <- range_bw[which(range_bw$gene %in% gene[g]),]
      plyr::ldply(seq_along(file_name),function(f){
        match_bw_rg <- tmp_g[which(tmp_g$fileName %in% file_name[f]),]$smax_value
        tmp <- junc_data[which(junc_data$fileName %in% file_name[f] &
                                 junc_data$gene %in% gene[g]),]
        tmp$dist <- scales::rescale(tmp$dist,to = c(0.5*match_bw_rg,0.9*match_bw_rg))
        
        return(tmp)
      }) -> tmp1
      return(tmp1)
    }) -> junc_data
    
    signal_layer_junction <- do.call(geom_arch,
                                     modifyList(
                                       list(data = junc_data,
                                            mapping = aes(x = start,xend = end,
                                                          color = fileName,
                                                          height = dist),
                                            linewidth = 0.5,
                                            show.legend = FALSE,
                                            guide = guide_legend(FALSE)),
                                       signal_layer_junction_params))
    
    
    
    signal_layer_junction_label <- do.call(geom_label,
                                           modifyList(
                                             list(data = junc_data,
                                                  mapping = aes(x = (start + end)/2,y = dist,
                                                                label = score),
                                                  label.size = NA),
                                             signal_layer_junction_label_params))
  }else{
    signal_layer_junction <- NULL
    signal_layer_junction_label <- NULL
  }
  
  
  
  # ==============================================================================
  # 13_main plot
  # ==============================================================================
  # color settings
  # sample_fill_col = sample_fill_col
  if (!is.null(bed.df)){
    peak_fill_col = RColorBrewer::brewer.pal(n = 12,name = "Paired")
    peak_fill_col = rep(c(peak_fill_col[seq(2,12,2)], peak_fill_col[seq(1,11,2)]), 4)[seq(length(unique(bed.df$sampleName)))]
  }
  
  # exon colors
  if(is.null(trans_fill_col)){
    trans_fill_col = RColorBrewer::brewer.pal(n = 9,name = "Paired")
  }else{
    trans_fill_col = trans_fill_col
  }
  
  useless_col = rep("white",1)
  
  
  # ============================================================
  # combine all layers
  pmain <-
    ggplot() +
    # bigwig layer
    ggnewscale::new_scale_fill() +
    signal_layer_geom_rect +
    scale_fill_manual(values = c(sample_fill_col,useless_col),
                      name = "") +
    # bed layer
    ggnewscale::new_scale_fill() +
    peaks_layer +
    scale_fill_manual(values = c(peak_fill_col,useless_col),name = "") +
    # junction layer
    ggnewscale::new_scale_colour() +
    signal_layer_junction +
    scale_color_manual(values = c(sample_fill_col,useless_col),
                       name = "") +
    signal_layer_junction_label +
    # other layers
    text_layer +
    seg_arrow_layer +
    # transcript layer
    trans_arrow_layer +
    ggnewscale::new_scale_fill() +
    trans_struct_layer +
    scale_fill_manual(values = trans_fill_col,name = "exon type") +
    theme_bw(base_size = base_size) +
    facet_layer +
    ggh4x::facetted_pos_scales(y = panel_range_layer,
                               x = panel_x_range_layer) +
    range_label_layer +
    panel_size_layer +
    theme(panel.grid = element_blank(),
          panel.spacing.y = unit(panel.spacing[2],"lines"),
          panel.spacing.x = unit(panel.spacing[1],"lines"),
          strip.text.y.left = element_text(angle = 0,face = "bold",hjust = 1),
          strip.text.x = element_text(face = "bold.italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          ggh4x.facet.nestline = element_line(colour = "black")) +
    xlab(xTitle) + ylab("")
  
  pmain
}




##
createSegment <- function(xPos = NULL,
                          yPos = NULL,
                          n_division = NULL,
                          rel_len = NULL,
                          abs_len = NULL){
  if(xPos[1] > xPos[2]){
    x = xPos[2];xend = xPos[1]
  }else{
    x = xPos[1];xend = xPos[2]
  }
  
  if(yPos[1] > yPos[2]){
    y = yPos[2];yend = yPos[1]
  }else{
    y = yPos[1];yend = yPos[2]
  }
  
  # get seqs
  if(y == yend){
    seg_len = abs(x - xend)
    if(!is.null(n_division)){
      seqs <- seq(from = x,to = xend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*seg_len
      }else{
        by = abs_len
      }
      seqs <- seq(from = x,to = xend,by = by)
      if(!(xend %in% seqs)) seqs <- c(seqs,xend)
    }
    res <- data.frame(id = 1:(length(seqs) - 1),
                      x = seqs[1:(length(seqs) - 1)],
                      xend = seqs[2:length(seqs)],
                      y = y,
                      yend = y)
  }else if(x == xend){
    seg_len = abs(y - yend)
    if(!is.null(n_division)){
      seqs <- seq(from = y,to = yend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*seg_len
      }else{
        by = abs_len
      }
      seqs <- seq(from = y,to = yend,by = by)
      if(!(yend %in% seqs)) seqs <- c(seqs,yend)
    }
    res <- data.frame(id = 1:(length(seqs) - 1),
                      x = x,
                      xend = x,
                      y = seqs[1:(length(seqs) - 1)],
                      yend = seqs[2:length(seqs)])
  }else{
    xseg_len = abs(xend - x)
    yseg_len = abs(yend - y)
    if(!is.null(n_division)){
      xseqs <- seq(from = x,to = xend,length.out = n_division)
      yseqs <- seq(from = y,to = yend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*xseg_len
      }else{
        by = abs_len
      }
      xseqs <- seq(from = x,to = xend,by = by)
      if(!(xend %in% xseqs)) xseqs <- c(xseqs,xend)
      yseqs <- seq(from = y,to = yend,length.out = length(xseqs))
    }
    res <- data.frame(id = 1:(length(xseqs) - 1),
                      x = xseqs[1:(length(xseqs) - 1)],
                      xend = xseqs[2:length(xseqs)],
                      y = yseqs[1:(length(yseqs) - 1)],
                      yend = yseqs[2:length(yseqs)])
  }
  return(res)
}

