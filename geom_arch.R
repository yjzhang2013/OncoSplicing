suppressMessages(library(GenomicRanges))


remove_args <- function(args, remove) {
  args[!names(args) %in% remove]
}

ggplot2_is_calculated_aes <- function(aesthetics){
  match <- "\\.\\.([a-zA-z._]+)\\.\\."
  stats <- rep(FALSE, length(aesthetics))
  grepl(match, sapply(aesthetics, deparse))
}

ggplot2_strip_dots <- function(aesthetics){
  match <- "\\.\\.([a-zA-z._]+)\\.\\."
  strings <- lapply(aesthetics, deparse)
  strings <- lapply(strings, gsub, pattern = match, replacement = "\\1")
  lapply(strings, function(x) parse(text = x)[[1]])
}


parseArgsForAes <- function(args){
  res <- if ("mapping" %in% names(args)) {
    args[["mapping"]]
  } else {
    Find(function(x) {
      class(x) == "uneval"
    }, args)
  }
  
  if (is.null(res))
    res <- structure(list(), class="uneval")
  
  idx <- ggplot2_is_calculated_aes(res)
  res[idx] <- ggplot2_strip_dots(res[idx])
  res
}

parseArgsForNonAes <- function(args){
  lst <- unlist(lapply(args, function(x) class(eval(x, parent.frame())) != "uneval"))
  args[lst]
}


setGeneric("geom_arch", function(data, ...) standardGeneric("geom_arch"))

setMethod("geom_arch", "data.frame", function(data, ..., 
                                              n = 25, max.height = 10){
  
  
  args <- list(...)
  args.aes <- parseArgsForAes(args)
  args.non <- parseArgsForNonAes(args)
  if("y" %in% names(args.aes))
    y.name <- quo_name(args.aes$y)
  else
    y.name <- NULL
  
  ## check required argument
  if(!all(c("x", "xend") %in% names(args.aes)))
    stop("x, xend, are requried in aes(), need to be passed into geom_arch()")
  startX <- rlang::eval_tidy(args.aes$x, data)
  endX <- rlang::eval_tidy(args.aes$xend, data)
  if("height" %in% names(args.aes)){
    if(!is.numeric(args.aes$height)){
      h <- rlang::eval_tidy(args.aes$height, data)
    }else{
      if(length(args.aes$height) == 1)
        h <- rep(args.aes$height, length(startX))
      else
        stop("unequal length of heights specified")
    }}else{
      h <- rep(max.height, length(startX))
    }
  if("y" %in% names(args.aes))
    y <- rlang::eval_tidy(args.aes$y, data)
  else
    y <- rep(0, length(startX))
  args.aes2 <- remove_args(args.aes, c("x", "y", "group", "hjust", "xend", "yend"))
  xx<-c()
  yy<-c()
  for(i in 1:n){
    ang<-i*pi/(2*n)
    xx[i]<-cos(ang)
    yy[i]<-sin(ang)
  }
  ##takes the quarter of the curve calculated, flips a copy over the y axis
  ##reduces time spent in for loop
  xx<-c(1,xx,rev(-xx),-1)
  yy<-c(0,yy,rev(yy), 0)
  ##SETS UP DATAFRAME TO KEEP TRACK OF ALL POINTS TO DRAW ALL ARCHES
  junc <- rep(seq_along(startX), each = length(xx))
  startX <- rep(startX, each = length(xx))
  endX <- rep(endX, each = length(xx))
  h <- rep(h, each = length(xx))
  y <- rep(y, each = length(xx))
  jump <- abs(endX - startX)
  jumpAdj <- if (length(jump)) max(jump) / max(abs(h)) else NA
  apoint <- data.frame(xx = xx * (abs(startX - endX) / 2) + (startX + endX) / 2,
                       yy = yy * h + y, junc,
                       s = ((abs(h) - jump / jumpAdj)) /
                         if (length(jump)) max(jump) else NA)
  data$junc <- seq_len(nrow(data))
  apoint <- merge(apoint, data, by = "junc")
  args.aes <- list(x = as.name("xx"),
                   y = as.name("yy"),
                   group = as.name("junc"))
  
  aesres <- do.call(aes, c(args.aes, args.aes2))
  aesres <- remove_args(aesres, "height")
  if(nrow(apoint)){
    args.non <- remove_args(args.non, "facets")
    reslst <- c(list(data = apoint), list(aesres),args.non)
    p <- do.call(geom_line, reslst)
    if("ylab" %in% names(args.non)){
      ylab <- quo_name(args.non$ylab)
    }else if(length(y.name)){
      ylab <- y.name
    }else{
      ylab <- ""
    }
    p <- list(p, ggplot2::ylab(ylab))
  }else{
    p <- NULL
  }
  p
})

## that means span the range of two end 
setMethod("geom_arch", "GRanges", function(data, ...,
                                           xlab, ylab, main,
                                           facets = NULL, rect.height = 0,
                                           n = 25, max.height = 10
){
  
  args <- list(...)  
  args$facets <- facets
  args.aes <- parseArgsForAes(args)
  args.non <- parseArgsForNonAes(args)
  args.non$max.height <- max.height
  facet <- build_facet(data, args)
  
  ## note rect.height = 0.4 is default cross ggbio
  ## need to make sure they are connected by two nearest point of rectangle
  df <- mold(data)
  if("height" %in% names(args.aes))
    signs <- sign(rlang::eval_tidy(args.aes$height, df))
  else
    signs <- 1
  args.aes$x <- substitute(start)
  args.aes$xend <- substitute(end)
  if("y" %in% names(args.aes)){
    y <- rlang::eval_tidy(args.aes$y, df)
    df[,quo_name(args.aes$y)] <- df[,quo_name(args.aes$y)] + rect.height * signs
  }else{
    df$.y <- rep(0, nrow(df)) + rect.height * signs
    args.aes$y <- substitute(.y)
  }
  if(nrow(df)){
    args.res <- c(list(data = df),
                  args.non,
                  list(do.call(aes, args.aes)))
    p <- do.call(geom_arch, args.res)
    p <- c(list(p) , list(facet))
  }else{
    p <- NULL
  }
  labels <- Labels(xlab, ylab, main, fallback = c(x = ""))
  p <- c(p, labels)
  p
})