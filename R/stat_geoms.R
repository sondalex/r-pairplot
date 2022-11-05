#' @export
StatPointRank <- ggproto("StatPointRank", StatIdentity,
  setup_data = function(data, params) {
    data$x = rank(data$x)
    data$y = rank(data$y)
    data
  }
)

#' @export
StatSmoothRank <- ggproto("StatSmoothRank", StatSmooth,
  setup_data = function(data, params) {
    data$x = rank(data$x)
    data$y = rank(data$y)
    data
  }
)

#' @export
geom_point_rank <- function(mapping = NULL, data = NULL, stat = "point_rank", position = "identity", 
    ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
    
    layer(data = data, mapping = mapping, stat = stat, geom = GeomPoint, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(na.rm = na.rm, ...))
}

#' @export
geom_smooth_rank <- function(mapping = NULL, data = NULL, stat = "smooth_rank", position = "identity", 
    ..., method = NULL, formula = NULL, se = TRUE, na.rm = FALSE, 
    orientation = NA, show.legend = NA, inherit.aes = TRUE) {
    params <- list(na.rm = na.rm, orientation = orientation, 
        se = se, ...)
    if (identical(stat, "smooth_rank")) {
        params$method <- method
        params$formula <- formula
    }
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSmooth, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = params)
}
