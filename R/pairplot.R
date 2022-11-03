#' @description Bin edges function based on 
#' https://github.com/numpy/numpy/blob/664a9f5cad40e958c73d2bc78f9bcd5bd5f818e6/numpy/lib/histograms.py
#' @param x A flat vector
#' @param range The range of bins.
#' @export
auto_bin <- function(x, range=NULL, na.rm=FALSE) {
  if (na.rm){
    x <- x[!(is.na(x))]
  }
  if (!length(range)){
    range <- c(min(x), max(x))
  } else {
    x = x[(x >= range[1]) & (x <= range[2]) ]
  }
  n = length(x)
  ptp <- (max(x) - min(x))
  sturges = ptp / (log2(n) + 1)
  fd = 2 * IQR(x)* n^(-1/3)
  auto <- min(sturges, fd)
  bincount <- ceiling((range[2] - range[1])/ auto)
  r <- seq(range[1], range[2], length.out=bincount + 1)
  return(r) 
}

#' @export
bin_width_auto <- function(x, na.rm=FALSE) {
  edges <- auto_bin(x, na.rm=na.rm)
  binwidth <- edges[2] - edges[1]
  return(binwidth)
}


#' Compute pairs of colnames based on the following architecture:
#' Suppose a dataframe with column ordered as follow : c("A", "B", "C", "D")
#' The pairs will be computed as follow
#' AA AB AC
#' BA BB BC
#' CA CB CC
#' @param data Dataframe
#' @return A list containing each pair. Would look like this 
#' [[1]]
#' [1] "A" "A"
#' 
#' [[2]]
#' [1] "A" "B"
#' 
#' [[3]]
#' [1] "A" "C"
#' 
#' [[4]]
#' [1] "B" "A"
#' 
#' [[5]]
#' [1] "B" "B"
#' 
#' [[6]]
#' [1] "B" "C"
#' 
#' [[7]]
#' [1] "C" "A"
#'
#' [[8]]
#' [1] "C" "B"
#' [[9]]
#' [1] "C" "C"
column_pairs <- function(data){
  columns <- colnames(data)
  cols_pairs <- vector(mode="list", length=length(columns)^2)
  k <- 1
  for (i in seq_along(columns)) {
    for (j in seq_along(columns)){
      pair <- c(columns[i], columns[j])
      cols_pairs[[k]] <- pair
      k <- k + 1
    }
  }
  return(cols_pairs)
}



add_modify_aes <- function(mapping, ...) {
  ggplot2:::rename_aes(modifyList(mapping, ...))  
}

# Operates on copy of list
#' @param xlim_func See [pairgrid()] parameter `common_xlim` (equivalent).
#' @param ylim_func A function to set y axis limits. If set in conjunction with `diag_share_ylim` set to TRUE.
#' Plots on the ylimit, will have their ylim set in common.
#' @param diag_share_ylim If set to TRUE. Y axis limits will also use `ylim_func`, except for diagonal only pairplot. In this case, 
#' will be normalized to unit length if the provided function implements unit scale normalization. See function
#' [pair_geom_histogram] for an example.
#' @param check.overlap See[ggplot2::guide_axis()]
#' @import ggplot2
#' @import rlang
mapplot <- function(data, mapping, pairs, grobs, indices, func, no_upper, no_lower, xlim_func=NULL, ylim_func=NULL, diag_share_ylim=TRUE, common_scale=NULL, check.overlap=TRUE, ...) {
  for (index in indices) {
    pair <- pairs[[index]]
    x <- pair[2]
    y <- pair[1]
    if (!length(func)) {
      # NULL was passed
      grobs[[index]] <- ggplot() + theme_void()
    } else {
        aes_map <- aes(x = !!rlang::sym(x), y=!!rlang::sym(y))
        if (length(mapping)) {
          mapping <- add_modify_aes(mapping, aes_map)
        } else {
          mapping <- aes_map
        }
        
        # If upper diagonal elements are NULL. Then 0-1 scale the diagonal.
        unit_y = ifelse(any(duplicated(pair)) & no_upper & no_lower, TRUE, FALSE)
        p <- func(data, mapping = mapping, unit_y=unit_y)
        if(length(xlim_func)) {
          xlim_pair <- xlim_func(data=data, x=x)
        } else{
          xlim_pair <- NULL
        }

        ylim_pair <- NULL
        if (length(ylim_func)) {
          if (diag_share_ylim & !unit_y) {
            ylim_pair=ylim_func(data=data, y=y)
          } else {
            # picky on selection: apply ylim to off diagonal only.
            if (x != y) {
              ylim_pair=ylim_func(data=data, y=y)
            }
          }
          p <- p + coord_cartesian(xlim = xlim_pair, ylim=ylim_pair)
        }
        # pass the identity guide.
        # ifelse on complicated object has uninented behaviour.
        guide <- if(check.overlap) guide_axis(check.overlap=check.overlap) else{waiver()}
        if (length(common_scale)) {
          # if values smaller than 1, display digits (for compatibility with scaled diagonal.)
          if (unit_y & any(duplicated(pair))) {
            p <- p + scale_y_continuous(labels=scales::label_number(accuracy=0.1), guide=guide)
          } else {
            p <- p + scale_y_continuous(labels=common_scale, guide=guide) 
          }
          p <- p + scale_x_continuous(labels=common_scale, guide=guide)
        } else {
          p <- p + scale_x_continuous(guide=guide) + scale_y_continuous(guide=guide)
        }
        
        grobs[[index]] <- p
      }
    }
  return(grobs)
}


#' @param data Data
#' @param x Column name for `x`. 
global_xlim <- function(data, x) {
  return(c(min(data[x], na.rm=T), max(data[x], na.rm=T)))
}

#' Description analog to [global_xlim()] but for y axis.
global_ylim <- function(data, y) {
  global_xlim(data=data, x=y)
}


#' Plot a pairplot.
#' @param size An integer representing the size of the square.
#' @param data A data frame. 
#' @param map_lower A function returning a ggplot2 object
#' @param map_diag Idem
#' @param map_upper Idem
#' @param only_x Boolean indicating whether the method passed to map_diag has aesteatics no y axis, such as histogram
#' @param common_xlim Accepts a function which returns a pair of limit. By default uses [global_xlim()]
#' @param diag_share_ylim See [mapplot()] for explanations.
#' @param common_scale A function for applying a common scale across all subplots.
#' @param repeat_text Boolean indicating whether the text (digits) should be set on all subplots. If FALSE, only the left and bottom contigous plots have text set.
#' @param ... Arguments passed to [mapplot()]
#'
#' @details Details 
#' 
#' ```{r}
#' 
#' library(ggplot2)
#' library(pairplot)
#' penguins_url <- "https://raw.githubusercontent.com/allisonhorst/palmerpenguins/main/inst/extdata/penguins.csv"
#' columns <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")
#' penguins <- read.csv(penguins_url)[columns]
#' penguins <- penguins[!apply(is.na(penguins), 1, any), ] # dropping NA

#' pairgrid(
#'           penguins,
#'           map_lower=pair_geom_smooth,
#'           map_diag=pair_geom_histogram,
#'           map_upper=pair_geom_point
#' )
#' ```
#' 
#' @export
#' @import ggplot2
#' @import patchwork
pairgrid <- function(data, mapping=NULL, map_lower, map_diag, map_upper, common_xlim=global_xlim, common_ylim=global_ylim, diag_share_ylim=TRUE, common_scale=scales::label_number(accuracy=1), repeat_labels=FALSE, repeat_text=FALSE, top=NULL, bottom=NULL, left=NULL, right=NULL, check.overlap=TRUE, ...) {
  size <- ncol(data)
  ncol <- size
  nrow <- ncol
  len <- ncol * nrow
  rows_col_pairs <- column_pairs(data)

  m <- matrix(seq(1, len), nrow, ncol, byrow = T)
  
  
  lower_index <- sort(m[lower.tri(m)])
  diag_index <- sort(diag(m))
  upper_index <- sort(m[upper.tri(m)])
   

  no_upper = is.null(map_upper)
  no_lower = is.null(map_lower)
  grobs <- vector(mode="list", length=len) # init
  
  grobs <- mapplot(
                   data, pairs=rows_col_pairs, grobs=grobs,
                   indices=lower_index, mapping = mapping,
                   func=map_lower, xlim_func = common_xlim,
                   ylim_func=common_ylim, diag_share_ylim=diag_share_ylim,
                   common_scale=common_scale, no_upper=no_upper, no_lower=no_lower,
                   check.overlap=check.overlap,
                   ...
                   )
  grobs <- mapplot(
                   data, pairs=rows_col_pairs, grobs=grobs,
                   indices=diag_index, mapping=mapping, func=map_diag,
                   xlim_func = common_xlim, ylim_func=common_ylim,
                   diag_share_ylim=diag_share_ylim,
                   common_scale=common_scale, no_upper=no_upper, no_lower=no_lower,
                   check.overlap=check.overlap,
                   ...
                   )
  grobs <- mapplot(
                   data, pairs=rows_col_pairs, grobs=grobs,
                   indices=upper_index, mapping=mapping, func=map_upper, 
                   xlim_func = common_xlim, ylim_func=common_ylim, 
                   diag_share_ylim=diag_share_ylim,
                   common_scale=common_scale,
                   no_upper=no_upper, no_lower=no_lower,
                   check.overlap=check.overlap,
                   ...
                   )
  
  vdelta <- ncol
  left_corner_index <- (vdelta * (nrow -1)) + 1
  left_contigous <- seq(1, left_corner_index, vdelta)
  bottom_contigous <- seq(left_corner_index, len)
  not_contigous <- setdiff(seq(1, len), union(bottom_contigous, left_contigous))
  
  for (index in seq_along(grobs)) {
    contigous_l <- index %in% left_contigous
    contigous_b <- index %in% bottom_contigous
    theme_i <- configure_theme(repeat_labels = repeat_labels, repeat_text = repeat_text, contigous_l = contigous_l, contigous_b = contigous_b)
    grobs[[index]] <- grobs[[index]] + theme_i
  }
  return(patchwork::wrap_plots(grobs, ncol=ncol, nrow=nrow, byrow=T))
}

#' Set theme programatically conditional on arguments.
#' @param repeat_label See [pairgrid()]
#' @param repeat_text See [pairgrid()]
#' @param contigous_l A boolean indicating whether subplot is contigous to the left side.
#' @param contigous_b A boolean indicating whether subplot is contigous to the bottom.
configure_theme <- function(repeat_labels, repeat_text, contigous_l, contigous_b) {
  general_theme <- theme_get()

  # I used switch in the past. It led to bug where NULL would be returned in some cases.
  axis.title.x <- if((!repeat_labels) & (!contigous_b)) element_blank() else{general_theme$axis.title.x}
  axis.title.y <- if((!repeat_labels) & (!contigous_l)) element_blank() else{general_theme$axis.title.y} 
  axis.text.x <- if((!repeat_text) & (!contigous_b)) element_blank() else{general_theme$axis.text.x}
  axis.text.y <- if((!repeat_text) & (!contigous_l)) element_blank() else{general_theme$axis.text.y}
  axis.ticks.x <- if((!repeat_text) & (!contigous_b)) element_blank() else{general_theme$axis.ticks.x}
  axis.ticks.y <- if((!repeat_text) & (!contigous_l)) element_blank() else{general_theme$axis.ticks.y}
  return(
    theme(
      axis.title.x = axis.title.x, axis.title.y=axis.title.y,
      axis.text.x = axis.text.x, axis.text.y = axis.text.y,
      axis.ticks.x = axis.ticks.x, axis.ticks.y=axis.ticks.y)
    )
}

test_configure_theme <- function() {
  # all combinations of argument
  
}

#' @export
pair_geom_smooth <- function(data, mapping, method="loess", formula=y~x, ...) {
  p <- ggplot(data = data, mapping = mapping) + 
      geom_point() + 
      geom_smooth(method=method, formula=formula, ...)
  return(p)
}

#' @export
pair_bar_plot <- function(data, mapping, stat="identity", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_bar(stat=stat, ...)
  return(p) 
}

#' @export
pair_geom_point <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) + 
      geom_point()
  return(p)
}


#' Min max histogram
#' @param unit_y Whether to scale y (min-max) to interval $[0, 1]$. If FALSE, min-max occurs but interval is set to $[min(y), max(y)]$, 
#' @export
pair_geom_histogram <- function(data, mapping, unit_y=FALSE, stat="bin", bins=10, binwidth=function(x) {bin_width_auto(x, na.rm=TRUE)}, ...) {
  bin_width_value <- binwidth(data[[1]])
  p <- ggplot(data = data, mapping = mapping)
  if (unit_y) {
    p <- p + geom_histogram(aes(y=..ncount..), stat=stat, bins=bins, binwidth=binwidth, ...)
  } else {
    # p <- p + geom_histogram(aes(y=after_stat(density / max(density, na.rm=T) * (max(x, na.rm=T) - min(x, na.rm=T) + min(x, na.rm=T))) , stat=stat, bins=bins, binwidth=binwidth, ...)
    p <- p + geom_histogram(aes(y=after_stat((max(x) - min(x)) * ((count) / (max(count) - min(count))) + min(x))) , stat=stat, bins=bins, binwidth=binwidth, ...)
  }
  return(p)
}

