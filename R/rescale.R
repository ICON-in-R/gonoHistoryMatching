
#' Take raw Latin Hypercube samples and rescale
#' 
#' @param ranges List of mean, sd for each parameter
#' @param vals Matrix of samples from Latin hypercube uniform distribution.
#'             Dimensions sample by parameter
#' @return Mean, upper and lower limit adjusted
#' 
rescale <- function(ranges, vals) {
  for (i in seq_along(ranges)) {
    rdiff <- ranges[[i]][2] - ranges[[i]][1]
    vals[, i] <- vals[, i]*rdiff + ranges[[i]][1] 
  }
  vals
}
