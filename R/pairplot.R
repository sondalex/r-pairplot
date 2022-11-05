#' @param i 
#' @param n
#' @return A named map: list(<subplot_index>=list(t=<int|NULL|>, b=<int|NULL|>, l=<int|NULL|>, r=<int|NULL|>)
neighbour_indices <- function(i, n, filter=NULL, upper_corner, lower_corner) {
  t <- i - n
  b <- i + n
  r <- i + 1
  l <- i - 1
  # replace by null values out of index.
  # And filter based on user filter
  if ((t < 1) || ((!("t" %in% filter)) & (i != n^2)) || (!upper_corner)){
    t = NULL
  }
  if ((b > n^2) || ((!("b" %in% filter)) & (i != 1)) || (!lower_corner)){
    b = NULL
  }
  if ((r > n^2) || ((!("r" %in% filter)) & (i != 1)) || (!upper_corner)) {
    r = NULL
  }
  if (l < 1 || (!("l" %in% filter) & (i != n^2)) || (!lower_corner)) {
    l = NULL
  }
  
  neighbours = vector(mode="list")
  neighbours[[as.character(i)]] = list(t=t, b=b, l=l, r=r)
  return(neighbours)
}

#' @return A named map `list(<subplot_index>=list(t=<ggdata|NULL>, b=<ggdata|NULL>, l=<ggdata|NULL>, r=<ggdata|NULL>)`
neighbour_data <- function(neighbour_pos, grobs){
  keys=ls(neighbour_pos)
  # values = unlist(neighbour_pos, use.names=F)
  newlist = list()
  for (key in keys){
  if (!is.null(neighbour_pos[[key]])){
    newlist[[key]] = ggplot_build(grobs[[neighbour_pos[[key]]]])
  } else {
    newlist[[key]] = NULL}
  }
  return(newlist)
}

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
#' In order to to construct a full pair matrix. This function has to be called
#' for each part (lower corner, upper corner, diagonal). The call on diagonal should be last
#' in order to have access to neighbour subplots metadata.
#' @param xlim_func See [pairgrid()] parameter `common_xlim` (equivalent).
#' @param ylim_func A function to set y axis limits. If set in conjunction with `diag_share_lim` set to TRUE,
#' plots on the ylimit, will have their ylim set in common.
#' @param diag_share_lim If set to TRUE. Y and X axis limits will also use `ylim_func`.
#' [pair_geom_histogram] for an example.
#' @param scale_diag_plot See [pairgrid()]
#' @param neighbours A map (named list) or `NULL`. 
#' Map of the form: `subplot_index=list(t=<int|NULL>, b=<int|NULL>, l=<int|NULL>, r=<int|NULL>)`
#' If not NULL, the function call
#' `func()` will include the data in argument `neighbour_ggdata` from the neighbour(s): `func(neighbour_ggdata=list(t=ggplot_build(grobs_index[[t_index]])), etc.)`.
#' There are no checks whether the selected grobs exist. Checks are reserved to function calling `mapplot`.
#' @param check.overlap See[ggplot2::guide_axis()]
#' @import ggplot2
#' @import rlang
mapplot <- function(data, mapping, pairs, grobs, indices, func, no_upper, no_lower, xlim_func=NULL, ylim_func=NULL, diag_share_lim=TRUE, common_scale=NULL, check.overlap=TRUE, neighbours=NULL, scale_diag_plot=TRUE, ...) {
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
        unit_y = ifelse(any(duplicated(pair)) & no_upper & no_lower & scale_diag_plot, TRUE, FALSE)
        if (length(neighbours)){
          current_neighbour_pos = neighbours[[as.character(index)]]
          near_data=neighbour_data(neighbour_pos=current_neighbour_pos, grobs = grobs)
          p <- func(data, mapping = mapping, unit_y=unit_y, neighbour_ggdata=near_data)
        } else {
          p <- func(data, mapping = mapping, unit_y=unit_y)
        }
        xlim_pair <- NULL
        if(length(xlim_func)){
          if (diag_share_lim) {
            xlim_pair <- xlim_func(data=data, x=x)
          } else {
            if (x != y) {
              xlim_pair=xlim_func(data=data, x=x)
            }
          }
        }

        ylim_pair <- NULL
        if (length(ylim_func)) {
          if (diag_share_lim) {
            ylim_pair=ylim_func(data=data, y=y)
          } else {
            # picky on selection: apply ylim to off diagonal only.
            if (x != y) {
              ylim_pair=ylim_func(data=data, y=y)
            }
          }
        }
        p <- p + coord_cartesian(xlim = xlim_pair, ylim=ylim_pair)
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
#' @param diag_share_lim See [mapplot()] for explanations.
#' @param common_scale A function for applying a common scale across all subplots.
#' @param repeat_text Boolean indicating whether the text (digits) should be set on all subplots. If FALSE, only the left and bottom contigous plots have text set.
#' @param text_on_diag Whether to apply text on diagonal. 
#' This argument applies only if one or more of the corners
#' have plots.
#' @param diag_neighbour A character vector c("t", "b", "l", "r"). If specified, 
#' the data from selected neighbour data is
#' passed to map_diag, functions. When specified the subplots 
#' on top-left and bottom-right corner will always return 
#' their neighbours data.
#' @param diag_share_ylim **Deprecated** Users should use diag_share_lim instead
#' @param scale_diag_plot Whether to unit scale diag only pairplot. If the function
#' passed to map_diag does not consider argument `unit_y` setting this option to TRUE might
#' lead to scale limits issues.
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
pairgrid <- function(data, mapping=NULL, map_lower, map_diag, map_upper, common_xlim=global_xlim, common_ylim=global_ylim, diag_share_lim=TRUE, common_scale=scales::label_number(accuracy=1), repeat_labels=FALSE, repeat_text=FALSE, check.overlap=TRUE, diag_neighbour=NULL, text_on_diag=T, diag_share_ylim=NULL, scale_diag_plot=FALSE, ...) {
  if (length(diag_share_ylim)) {
    stop("diag_share_ylim is deprecated and has been replaced by diag_share_lim.")
  }
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
  
  if (length(diag_neighbour)) {
    if (any(!(diag_neighbour %in% c("t", "l", "b", "r")))) {
      stop(
        paste0("Characters must be in c('t', 'b', 'l', 'r')): ")
      )
    }
    # data should be accessible via 
    # neighbour_ggdata$t, neighbour_ggdata$b, 
    # neighbour_ggdata$l, neighbour_ggdata$r
    # accessible via t, b, l, r.
    
    map_neighbours = sapply(diag_index, function(i){
      neighbour_indices(i, n=size, filter=diag_neighbour, lower_corner=!is.null(map_lower), upper_corner=!is.null(map_upper))
      }
    )
  } else {
    map_neighbours = NULL
  }
  

  grobs <- mapplot(
                   data, pairs=rows_col_pairs, grobs=grobs,
                   indices=lower_index, mapping = mapping,
                   func=map_lower, xlim_func = common_xlim,
                   ylim_func=common_ylim, diag_share_lim=diag_share_lim,
                   common_scale=common_scale, no_upper=no_upper, no_lower=no_lower,
                   check.overlap=check.overlap, neighbour=NULL, scale_diag_plot=scale_diag_plot,
                   ...
                   )
  grobs <- mapplot(
                   data, pairs=rows_col_pairs, grobs=grobs,
                   indices=upper_index, mapping=mapping, func=map_upper, 
                   xlim_func = common_xlim, ylim_func=common_ylim, 
                   diag_share_lim=diag_share_lim,
                   common_scale=common_scale,
                   no_upper=no_upper, no_lower=no_lower,
                   check.overlap=check.overlap, neighbour=NULL,
                   scale_diag_plot=scale_diag_plot,
                   ...
                   )
  grobs <- mapplot(
                   data, pairs=rows_col_pairs, grobs=grobs,
                   indices=diag_index, mapping=mapping, func=map_diag,
                   xlim_func = common_xlim, ylim_func=common_ylim,
                   diag_share_lim=diag_share_lim,
                   common_scale=common_scale, no_upper=no_upper, no_lower=no_lower,
                   check.overlap=check.overlap, neighbour=map_neighbours,
                   scale_diag_plot=scale_diag_plot,
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
    on_diag <- ifelse(index %in% diag_index, TRUE, FALSE)
    theme_i <- configure_theme(repeat_labels = repeat_labels, repeat_text = repeat_text, contigous_l = contigous_l, contigous_b = contigous_b, on_diag=on_diag, text_on_diag=text_on_diag)
    grobs[[index]] <- grobs[[index]] + theme_i
  }
  return(patchwork::wrap_plots(grobs, ncol=ncol, nrow=nrow, byrow=T))
}

#' Set theme programatically conditional on arguments.
#' @param repeat_label See [pairgrid()]
#' @param repeat_text See [pairgrid()]
#' @param text_on_diag See [pairgrid()]
#' @param contigous_l A boolean indicating whether subplot is contigous to the left side.
#' @param contigous_b A boolean indicating whether subplot is contigous to the bottom.
#' @param on_diag A boolean indicating whether the subplot is on diagonal.
configure_theme <- function(repeat_labels, repeat_text, text_on_diag, contigous_l, contigous_b, on_diag) {
  general_theme <- theme_get()
  diag = text_on_diag & on_diag # whether effective on diag
  # I used switch in the past. It led to bug where NULL would be returned in some cases.
  axis.title.x <- if((!repeat_labels) & (!contigous_b)) element_blank() else{general_theme$axis.title.x}
  axis.title.y <- if((!repeat_labels) & (!contigous_l)) element_blank() else{general_theme$axis.title.y} 
  axis.text.x <- if((!repeat_text) & (!contigous_b) & (!diag)) element_blank() else{general_theme$axis.text.x}
  axis.text.y <- if((!repeat_text) & (!contigous_l) & (!diag)) element_blank() else{general_theme$axis.text.y}
  axis.ticks.x <- if((!repeat_text) & (!contigous_b) & (!diag)) element_blank() else{general_theme$axis.ticks.x}
  axis.ticks.y <- if((!repeat_text) & (!contigous_l) & (!diag)) element_blank() else{general_theme$axis.ticks.y}
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

#' @x A vector
min_max_scale <- function(x, a, b) {
  (b - a) * ((x - min(x)) / (max(x) - min(x))) + a
}


#' Min max histogram
#' @param unit_y Whether to scale y (min-max) to interval [0, 1].
#' @param min_max A function for transforming the data.
#' Those values are used as boundaries. If only x or y is provided. Scaling is done on one dimension only.
#' for min max scaling.
#' @export
pair_geom_histogram <- function(data, mapping, neighbour_ggdata=NULL, unit_y=FALSE, stat="bin", bins=10, binwidth=function(x) {bin_width_auto(x, na.rm=TRUE)}, ...) {
  bin_width_value <- binwidth(data[[1]])
  x = rlang::as_name(mapping$x)
  if (unit_y) {
    p <- ggplot(data = data, mapping = mapping) + geom_histogram(aes(y=..ncount..), stat=stat, bins=bins, binwidth=binwidth, ...)
  } else {
    # p <- p + geom_histogram(aes(y=after_stat(density / max(density, na.rm=T) * (max(x, na.rm=T) - min(x, na.rm=T) + min(x, na.rm=T))) , stat=stat, bins=bins, binwidth=binwidth, ...)
    if (length(neighbour_ggdata)){
        for (key in c("t", "b")) {
          xdata = neighbour_ggdata[[key]]
          if (length(xdata)){
            break
          }
        }
        for (key in c("l", "r")) {
          ydata = neighbour_ggdata[[key]]
          if (length(ydata)){
            break
          }
        }
        
        # Two corner case
        # 1)
        # 1 NA NA 
        # 1  1 NA
        # 1  1  1

        # --> TL has ydata NULL and xdata existing 
        # --> BR has ydata existing, but xdata NULL

        # 2)
        # 1  1  1
        # NA 1  1
        # NA NA 1
        
        # --> TL has ydata existing, but xdata not existing.
        # --> BR has ydata NULL but xdata existing.
        
        # The bounds for axis which does not share a neighbour
        # are set to the minimum of the data.

        # Combination of upper_corner and lower corner have to be considered
        # for the top left and bottom right graphic.
        # get first non-null for y and x axis.
        
        xaxis_bounds <- if(length(xdata)) xdata$layout$panel_scales_x[[1]]$range$range else {c(min(data[x]), max(data[x]))}
        yaxis_bounds <- if(length(ydata)) ydata$layout$panel_scales_y[[1]]$range$range else {c(min(data[x]), max(data[x]))}
        
        lower_bound_y <- yaxis_bounds[1]
        upper_bound_y <- yaxis_bounds[2]
        lower_bound_x <- xaxis_bounds[1]
        upper_bound_x <- xaxis_bounds[2]
        
        t_data = min_max_scale(data[x], lower_bound_x, upper_bound_x)
        p <- ggplot(data = t_data, mapping = mapping)
        # min max-scale a first time. Such that the x-axis is scaled.
    } else {
      p <- ggplot(data = data, mapping = mapping)
      upper_bound_y <- max(data[x])
      lower_bound_y <- min(data[x]) 
    }
    p <- p + geom_histogram(aes(y=after_stat((upper_bound_y - lower_bound_y) * ((count - min(count)) / (max(count) - min(count))) + lower_bound_y)), stat=stat, bins=bins, binwidth=binwidth, ...)
  }

  return(p)
}
